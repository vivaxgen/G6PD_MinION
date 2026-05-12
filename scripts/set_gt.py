import pathlib
from cyvcf2 import VCF, Writer
import numpy as np
import argparse
from ngs_pipeline import cout, cerr, cexit

ap = argparse.ArgumentParser(description="Set GT fields in VCF file")
ap.add_argument("--infile", help="input VCF file")
ap.add_argument("--outfile", help="output VCF file")
ap.add_argument("--minimum_depth", help="minimum depth to call a genotype", type=int, default=-1)
ap.add_argument("--minimum_minor_depth", help="minimum minor allele depth to call a heterozygous genotype", type=int, default=-1)
ap.add_argument("--minimum_minor_ratio", help="minimum minor allele ratio to call a heterozygous genotype", type=float, default=-1)
ap.add_argument("--set_het_to_ref", help="set heterozygous genotypes to reference", action="store_true")
ap.add_argument("--set_het_to_alt", help="set heterozygous genotypes to alternate", action="store_true")
ap.add_argument("--set_het_to_missing", help="set heterozygous genotypes to missing", action="store_true")
ap.add_argument("--set_missing_to_het", help="set missing genotypes to heterozygous", action="store_true")
ap.add_argument("--set_missing_to_ref", help="set missing genotypes to reference", action="store_true")
ap.add_argument("--set_missing_to_alt", help="set missing genotypes to alternate", action="store_true")
ap.add_argument("--set_alt2_to_ref", help="set genotypes with alternate allele index > 1 to reference", action="store_true")
ap.add_argument("--set_alt2_to_alt1", help="set genotypes with alternate allele index > 1 to alternate allele index 1", action="store_true")
ap.add_argument("--set_alt2_to_missing", help="set genotypes with alternate allele index > 1 to missing", action="store_true")
ap.add_argument("--set_id", help="set ID field to CHROM:POS", action="store_true")
ap.add_argument("--headers", help="additional header line to add to the output VCF file", action="append", default=[])
ap.add_argument("--threads", help="number of threads to use for reading VCF file", type=int, default=1)


def set_genotypes(variant, indexes, value):
    for idx in indexes:
        variant.genotypes[idx] = value


def set_alt2_gt(variant, allele):
    # find positions that have allele index > 1
    # for idx in range(len(variant.genotypes)):
    #    if variant.genotypes[idx][0]
    for idx in range(len(variant.genotypes)):
        if variant.genotypes[idx][0] > 1:
            variant.genotypes[idx] = [
                allele,
                variant.genotypes[idx][1],
                variant.genotypes[idx][2],
            ]
        if variant.genotypes[idx][1] > 1:
            variant.genotypes[idx] = [
                variant.genotypes[idx][0],
                allele,
                variant.genotypes[idx][2],
            ]

def set_GT(
    infile: str | pathlib.Path,
    outfile: str | pathlib.Path,
    *,
    minimum_depth: int = -1,
    minimum_minor_depth: int = -1,
    minimum_minor_ratio: float = -1,
    set_het_to_ref: bool = False,
    set_het_to_alt: bool = False,
    set_het_to_missing: bool = False,
    set_missing_to_het: bool = False,
    set_missing_to_ref: bool = False,
    set_missing_to_alt: bool = False,
    set_alt2_to_ref: bool = False,
    set_alt2_to_alt1: bool = False,
    set_alt2_to_missing: bool = False,
    set_id: bool = False,
    headers: list = [],
    threads: int = 1,
) -> None:

    genotype_HOM_REF = [0, 0, False]
    genotype_HET = [0, 1, False]
    genotype_MISSING = [-1, -1, False]
    genotype_HOM_ALT = [1, 1, False]

    vcf = VCF(infile, threads=threads)
    w = Writer(outfile, vcf)

    for header in headers:
        w.add_to_header(header)

    logs = []

    USE_GT = True
    if minimum_minor_depth >= 0 or minimum_minor_ratio >= 0:
        if minimum_depth < 0:
            raise ValueError(
                "when setting for minimum_minor_depth or minimum_minor_ratio, "
                "minimum_depth also needs to be set"
            )
        USE_GT = False


    for v in vcf:

        AD = v.format("AD")
        if AD is None: #clair3 sometimes do ./. with only GT for targeted panel
            AD = np.zeros((len(v.genotypes), 2), dtype=int)
        AD = AD.clip(min=0) # handle ".,." returned as negative number 

        if not USE_GT:
            if len(v.ALT) == 0:
                # in case no alternate base
                major_alleles = minor_alleles = minor_depths = np.zeros(AD.shape[0], dtype=int)
                major_depths = AD[:, 0]
            else:
                # we have biallelic alleles or multiple alternate alleles, but only care about 2 alleles (major and minor)
                gt_argsorted = np.argsort(-AD, stable = True, axis=1) #argsort is stable, when equal
                major_alleles = gt_argsorted[:, 0]
                minor_alleles = gt_argsorted[:, 1]

                major_depths = AD[np.arange(AD.shape[0]), major_alleles]
                minor_depths = AD[np.arange(AD.shape[0]), minor_alleles]

            # up to this point, there should be:
            # major_alleles
            # minor_alleles (after adjustment)
            # minor_depths (after adjustment)

            minor_ratios = minor_depths / (major_depths + minor_depths)

            non_hets = (minor_depths < minimum_minor_depth) | (
                minor_ratios < minimum_minor_ratio
            )

            # set minor_alleles to major_alleles if minor_depths == 0;
            null_minor_indexes = minor_depths == 0
            minor_alleles[null_minor_indexes] = major_alleles[null_minor_indexes]

            # for all that are non-hets, set both alleles to major allele
            for idx in non_hets.nonzero()[0]:
                allele = major_alleles[idx]
                v.genotypes[idx] = [allele, allele, v.genotypes[idx][-1]]

            # hets are those that are not non-hets
            hets = (~non_hets).nonzero()[0]
            non_hets = non_hets.nonzero()[0]
        
        else:
            # get non-hets from genotypes / GT fields
            hets = [i for i, gt in enumerate(v.genotypes) if (len(set(gt[:-1])) == 2 and not -1 in gt)]
            non_hets = [i for i, gt in enumerate(v.genotypes) if (len(set(gt[:-1])) == 1 and not -1 in gt)]

            gts = np.array(v.genotypes)[:, :2]
            gt0_AD = AD[np.arange(AD.shape[0]), gts[:, 0]]
            gt1_AD = AD[np.arange(AD.shape[0]), gts[:, 1]]
            gt1_AD[non_hets] = 0

            major_alleles = np.where(gt0_AD >= gt1_AD, gts[:, 0], gts[:, 1])
            minor_alleles = np.where(gt0_AD < gt1_AD, gts[:, 0], gts[:, 1])
            major_depths = np.where(gt0_AD >= gt1_AD, gt0_AD, gt1_AD)
            minor_depths = np.where(gt0_AD < gt1_AD, gt0_AD, gt1_AD)


        gt_depths = major_depths + minor_depths


        # -- handling hets -- #

        # for all hets, set their alleles
        if set_het_to_ref:
            set_genotypes(v, hets, genotype_HOM_REF)
        elif set_het_to_alt:
            set_genotypes(v, hets, genotype_HOM_ALT)
        elif set_het_to_missing:
            set_genotypes(v, hets, genotype_MISSING)
        elif minimum_minor_depth >= 0 or minimum_minor_ratio >= 0:
            for idx in hets:
                # sort alleles to give 0/1, 0/2 or 1/2 instead of 1/0, 2/0 or 2/1
                maj_allele = major_alleles[idx]
                min_allele = minor_alleles[idx]
                alleles = [*sorted([min_allele, maj_allele]), False]
                v.genotypes[idx] = alleles

        # -- handling missing -- #

        if minimum_depth >= 0:
            missings = (gt_depths < minimum_depth).nonzero()[0]
        else:
            missings = [i for i, gt in enumerate(v.genotypes) if -1 in gt]

        if set_missing_to_ref:
            set_genotypes(v, missings, genotype_HOM_REF)
        elif set_missing_to_alt:
            # WARNING: this is correct only for biallelic variants
            set_genotypes(v, missings, genotype_HOM_ALT)
        elif set_missing_to_het:
            # WARNING: this is correct only for biallelic variants
            set_genotypes(v, missings, genotype_HET)
        elif minimum_depth >= 0:
            set_genotypes(v, missings, genotype_MISSING)

        # handling other alternate alleles

        if set_alt2_to_ref:
            set_alt2_gt(v, allele=0)
        elif set_alt2_to_alt1:
            set_alt2_gt(v, allele=1)
        elif set_alt2_to_missing:
            set_alt2_gt(v, allele=-1)
            # XXX: need to check if alleles = ./? or ?/., then set to ./.

        # reset genotypes
        v.genotypes = v.genotypes

        if set_id:
            v.ID = f"{v.CHROM}:{v.POS}"

        w.write_record(v)

    vcf.close()
    w.close()

    cerr(f"Processsed VCF file written to {outfile}")

if __name__ == "__main__":
    args = ap.parse_args()
    set_GT(**vars(args))