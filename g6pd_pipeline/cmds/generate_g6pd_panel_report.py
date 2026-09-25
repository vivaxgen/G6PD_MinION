import argparse
from ngs_pipeline import cerr, cexit

def init_argparser():
    p = argparse.ArgumentParser(description='generate variant reports from VCF')
    p.add_argument('--infofile', required=True,
                   help='a TSV file containing variant information')
    p.add_argument('-o', '--outfile', required=True,
                   help='output filename, eg. my-report.tsv')
    p.add_argument('--mindepth', type=int, default=20,
                   help='minimum depth for calling variants')
    p.add_argument('--min_var_qual', type=float, default=10,
                   help='minimum quality for variants')
    p.add_argument('infile',
                   help='input file in vcf.gz format')
    p.add_argument('--clair3_gvcf', action='store_true',
                   help='input file is in Clair3 gvcf format')
    p.add_argument('--no_flag_failed_variant', action='store_false', default=True,
                   help='Output all variant, including those with depth < mindepth marked with (*) and qual < minqual marked with (^)')
    p.add_argument("--multiple_missenses", type=str, default=None,
                   help="TSV containing the pseudo-phased multiple missense variants")
    return p

def gts_to_gt_changes(gts, ref, alts):
    changes = []
    for gt in gts:
        if gt == -1:
            changes.append("?")
        elif gt == 0:
            changes.append(ref + "=")
        else:
            changes.append(ref + ">" + alts[gt-1])
    return set(changes)

def generate_variant_report(args):
    import pandas as pd
    import sys
    from cyvcf2 import VCF
    # import IPython

    cerr(f'[Reading variant information from {args.infofile}]')

    info_df = pd.read_table(args.infofile)
    info_df.set_index("ID", inplace=True)
    info_df = info_df.sort_index().sort_values(by=["CHROM", "POS"])
    grouped_names = (
        info_df.groupby("ID", sort=False)
        .apply(
            lambda group: (
                group["Name"].iloc[0]
                + " ("
                + " + ".join(
                    group["CHROM"].astype(str) + ":" + group["POS"].astype(str)
                )
                + ")"
            ),
            include_groups=False,
        )
    )
    info_df["var_name"] = info_df.index.map(grouped_names)
    info_df = info_df.reset_index()
    cerr(f'[Generating variant report from {args.infile}]')

    variants_df_list = []
    vcf = VCF(args.infile, gts012=True)
    sample = vcf.samples[0]
    is_clair3_gvcf = args.clair3_gvcf
    for v in vcf:
        try:
            if is_clair3_gvcf:
                # if GT != 0/0, assert that the ALT is not <NON_REF>
                ALT = max(v.genotypes[0][0:2])
                # Sanity check, if failed, something must have gone wrong
                if ALT > 0:
                    assert([v.REF, *v.ALT][ALT] != "<NON_REF>")

            if 'DP' in v.FORMAT:
                allele_depth = v.format('DP')[0][0]
            elif 'MIN_DP' in v.FORMAT:
                allele_depth = v.format('MIN_DP')[0][0]
            else:
                try:
                    allele_depth = v.INFO['DP']
                except KeyError:
                    allele_depth = 0
            
            chrom = v.CHROM
            pos = v.POS
            ref = v.REF
            alt = v.ALT
            depth = max(0, allele_depth)
            qual = max(0, v.QUAL) if v.QUAL is not None else 0
            gt = v.genotypes[0][:-1]
            ad = v.format('AD')[0].clip(0) if 'AD' in v.FORMAT else [0, 0]

            variants_df_list.append(pd.DataFrame([[chrom, pos, ref, alt, depth, qual, gt, ad]]))

        except KeyError:
            cerr(f'[WARNING] is_clair_gvcf:{is_clair3_gvcf} - {v.CHROM}:{v.POS} '
                 f'REF: {v.REF}, ALT: {v.ALT} is not found in the infofile, skipping')

    if len(variants_df_list) > 0:
        variants_df = pd.concat(variants_df_list, ignore_index=True)
        variants_df.columns = ["CHROM", "POS", "REF", "ALT", "DP", "QUAL", "GT", "AD"]
    else:
        variants_df = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT", "DP", "QUAL", "GT", "AD"])

    # handle freebayes duplicates
    # duplicated entry with all blank values, remove them
    dup_entries = variants_df[["CHROM", "POS"]].value_counts().to_frame("count").query("count > 1").reset_index()
    idx_to_drop = []
    for _, row in dup_entries.iterrows():
        chrom = row["CHROM"]
        pos = row["POS"]
        dup_rows = variants_df.query("CHROM == @chrom and POS == @pos").copy()
        dup_rows["alt_s"] = dup_rows["ALT"].apply(lambda x: ",".join(set(sorted(x))))
        if dup_rows["REF"].nunique() == 1 and dup_rows["alt_s"].nunique() == 1:
            # get row with DP otherwise, keep the first row
            if dup_rows["DP"].sum() > 0:
                idx_to_drop.extend(dup_rows.query("DP == 0").index.tolist())
            else:
                idx_to_drop.extend(dup_rows.index.tolist()[1:])

    variants_df = variants_df.drop(idx_to_drop).reset_index(drop=True)
    # fill in missing for variants_df
    missing_from_vcf = info_df[~info_df.set_index(["CHROM", "POS"]).index.isin(variants_df.set_index(["CHROM", "POS"]).index)]
    fill_in_df_list = []
    for _, row in missing_from_vcf.iterrows():
        fill_in_df_list.append(
            pd.DataFrame({
                "CHROM": row["CHROM"],
                "POS": row["POS"],
                "REF": row["Change"].split(">")[0],
                "ALT": [[]],
                "DP": 0,
                "QUAL": 0,
                "GT": [[-1, -1]],
                "AD": [[0]],
            }, index=[0])
        )
    variants_df = pd.concat([variants_df, *fill_in_df_list], ignore_index=True)

    result_df = info_df.copy()
    known_alleles = result_df.groupby(["CHROM", "POS"])["Change"].apply(set).reset_index()
    known_alleles["known_alts"] = known_alleles["Change"].apply(lambda x: [a.split(">")[1] for a in x])
    result_df = result_df.merge(known_alleles[["CHROM", "POS", "known_alts"]], on=["CHROM", "POS"], how="left")

    # fill in from variants_df
    result_df = result_df.merge(variants_df, left_on=["CHROM", "POS"], right_on=["CHROM", "POS"], how="outer")
    result_df["gt_change"] = result_df.apply(lambda row: gts_to_gt_changes(row['GT'], row['REF'], row['ALT']), axis=1)

    def calculate_unknown_allele_depth_ratio(row):
        alts_idx = [i for i, alt in enumerate(row["ALT"]) if alt in row["known_alts"]]
        known_allele_depth = 0
        if row["DP"] > 0:
            if row["AD"] is not None:
                if len(row["AD"]) > 0:
                    known_allele_depth = row["AD"][0] + sum([row["AD"][i+1] for i in alts_idx])
            return 1 - known_allele_depth / row["DP"]
        else:
            return 1

    result_df["unknown_base"] = result_df.apply(lambda row: calculate_unknown_allele_depth_ratio(row), axis = 1)

    def fill_found_variants_single(row):
        assert len(row["gt_change"]) <= 2, f"Unexpected gt_change value: {row['gt_change']}"
        if row["gt_change"] == {'?'}:
            variant = "?"
        # home/hemi
        elif len(row["gt_change"]) == 1:
            if row["Change"] in row["gt_change"]: # interested GT change is found
                variant = "+"
            else:
                variant = "-"
        # het
        else:
            if row["Change"] in row["gt_change"]: # interested GT change is found
                variant = "-/+"
            else:
                variant = "-"
        
        if not args.no_flag_failed_variant:
            return variant
        else:
            failed_depth = ""
            failed_qual = ""
            if variant != "?":
                # het threshold and depth for alleles handled upstream
                if row["DP"] < args.mindepth:
                    failed_depth = "d"
                if row["QUAL"] < args.min_var_qual:
                    failed_qual = "q"
            return f"{variant}{failed_depth}{failed_qual}"    
        
    result_df.loc[:, sample] = result_df.apply(lambda row: fill_found_variants_single(row), axis=1)
    n_failed = result_df.loc[result_df["DP"] < args.mindepth, ["CHROM", "POS"]].drop_duplicates().shape[0]

    sample_is_all_hemi_homo = all(result_df[sample].isin(["+", "-", "+d", "-d", "+q", "-q", "+dq", "-dq"]))

    multiple_missenses_df = pd.read_table(args.multiple_missenses) if args.multiple_missenses else pd.DataFrame()
    if not multiple_missenses_df.empty:
        multiple_missenses_df["ratio"] = multiple_missenses_df.groupby("marker")["count"].transform(lambda x: x / x.sum())
        non_unknown_counts = (
            multiple_missenses_df.loc[
                ~multiple_missenses_df["haplotype"].str.contains("?", regex=False),
                ["marker", "count"],
            ]
            .groupby("marker")["count"]
            .sum()
        )
        multiple_missenses_df["ratio_non_unknown"] = multiple_missenses_df["count"] / multiple_missenses_df["marker"].map(non_unknown_counts)

    multiple_missenses_variants = result_df["ID"].value_counts().loc[lambda x: x > 1].index.tolist()


    def confirm_phase_result(rows, check_phased):
        marker = rows["ID"].values[0]
        default = ("u", 0, 0)
        if not multiple_missenses_df.empty:
            if marker in multiple_missenses_df["marker"].values:
                phasing_result = multiple_missenses_df.query("marker == @marker")
                if not phasing_result.empty:
                    found = check_phased in phasing_result["haplotype"].values
                    if found:
                        result = phasing_result.query("haplotype == @check_phased")
                        return ("", result["ratio"].values[0], result["ratio_non_unknown"].values[0])
        return default


    def get_phase_result(rows, positive):
        marker = rows["ID"].values[0]
        if not multiple_missenses_df.empty:
            if marker in multiple_missenses_df["marker"].values:
                phasing_result = multiple_missenses_df.query("marker == @marker and not haplotype.str.contains('?', regex=False)")
                if not phasing_result.empty:
                    trustable = ""
                    phasing_result = phasing_result.sort_values(by="ratio_non_unknown", ascending=False)
                    if phasing_result["ratio"].sum() < 0.6:
                        trustable = "#"
                    if phasing_result.shape[0] > 2:
                        trustable = "#"
                    if positive in phasing_result.head(2)["haplotype"].values:
                        return "-/+", trustable
                    else:
                        return "-", trustable
        return "-/+", "u"

    # final result = phased result + " ; " + single result
    # ? if any of the multiple missense variants is ?
    # u if unphaseable
    # # if phased result is + but the ratio of phased variants is < 0.8 (based on vcf)
    # - if none of the multiple missense variants is + or -/+
    # + if all of the multiple missense variants were + and they are phased together
    # or 
    # -/+ if both of the multiple missense variants were -/+ and they are phased together
    # if any of the multiple missense variants were -/+ and the others are + and they are phased together, then -/+
    for multiple_missenses_var in multiple_missenses_variants:
        interested_rows = result_df.loc[result_df["ID"] == multiple_missenses_var, :].copy()
        interested_rows = interested_rows.sort_values(by=["CHROM", "POS"])
        required_phase_result = "".join(interested_rows["Change"].str.split(">").str[-1])

        if any(interested_rows[sample].isin(["?", "?d", "?q", "?dq"])):
            phase_result = "?"
            phase_result_confirmation = ""
        elif any(interested_rows[sample].isin(["-", "-d", "-q", "-dq"])):
            phase_result = "-"
            phase_result_confirmation = ""
        elif all(interested_rows[sample].isin(["+", "+d", "+q", "+dq"])):
            phase_result = "+"
            phase_result_confirmation, phase_ratio, phase_ratio_non_unknown = confirm_phase_result(interested_rows, check_phased = required_phase_result)
            if phase_ratio_non_unknown < 0.7 and phase_result_confirmation == "":
                phase_result_confirmation = "#"
        else: ### Only handling het cases here
            phase_result, phase_result_confirmation = get_phase_result(interested_rows, required_phase_result)

        final_phase_result = f"{phase_result}{phase_result_confirmation}"
        result_df.loc[result_df["ID"] == multiple_missenses_var, sample] =  phase_result + " ; " + result_df.loc[result_df["ID"] == multiple_missenses_var, sample]

    interested_columns = [
        "CHROM", "POS", "Change", "Name", "Amplicon", "var_name", "REF", "ALT", "DP", "QUAL", "GT", "AD", "gt_change", "unknown_base", sample
    ]
    full_details_tsv_name = args.outfile.replace(".tsv", ".full_details.tsv")
    result_df[interested_columns].to_csv(full_details_tsv_name, index=False, sep='\t')

    pos_hom_var = []
    pos_het_var = []
    neg_var = []
    pos_failed_hom_var = []
    pos_failed_het_var = []
    neg_failed_var = []
    pos_hom_multiple_missense_unphaseable = []
    pos_hom_multiple_missense_lowconf = []
    pos_het_multiple_missense_unphaseable = []
    pos_het_multiple_missense_lowconf = []
    unknown_var = []

    for var in result_df["var_name"].unique():
        var_rows = result_df.loc[result_df["var_name"] == var, :].sort_values(by=["CHROM", "POS"])
        var_result = var_rows[sample].values[0]
        if var_rows.shape[0] == 1:
            var_result = var_rows[sample].values[0]
            match var_result:
                case "+":
                    pos_hom_var.append(var)
                case "-" | "-q" :
                    neg_var.append(var)
                case "-/+":
                    pos_het_var.append(var)
                case "?":
                    unknown_var.append(var)
                case "+d" | "+q" | "+dq":
                    pos_failed_hom_var.append(var)
                case "-/+d" | "-/+q" | "-/+dq":
                    pos_failed_het_var.append(var)
                case "-d" | "-dq":
                    neg_failed_var.append(var)
                case "?":
                    unknown_var.append(var)
                case _:
                    raise ValueError(f"Unexpected variant result: {var_result} for {var}")
        else:
            phased_result = var_rows[sample].str.split(" ; ").str[0].values[0]
            ind_result = ", ".join(var_rows[sample].str.split(" ; ").str[1].values)
            match phased_result:
                case "+":
                    pos_hom_var.append(var)
                case "-":
                    neg_var.append(var)
                case "-/+":
                    pos_het_var.append(var)
                case "?":
                    unknown_var.append(var)
                case "+u":
                    pos_hom_multiple_missense_unphaseable.append(var)
                case "+#":
                    pos_hom_multiple_missense_lowconf.append(var)
                case "-/+u":
                    pos_het_multiple_missense_unphaseable.append(var)
                case "-/+#":
                    pos_het_multiple_missense_lowconf.append(var)
                case _:
                    raise ValueError(f"Unexpected phased result: {phased_result} for {var}")

    report_df = pd.DataFrame({
        "sample": sample,
        "pos_hom_var": "\n".join(pos_hom_var),
        "pos_het_var": "\n".join(pos_het_var),
        "neg_var": "\n".join(neg_var),
        "pos_failed_het_var": "\n".join(pos_failed_het_var),
        "pos_failed_hom_var": "\n".join(pos_failed_hom_var),
        "neg_failed_var": "\n".join(neg_failed_var),
        "pos_hom_multiple_missense_unphaseable": "\n".join(pos_hom_multiple_missense_unphaseable),
        "pos_hom_multiple_missense_lowconf": "\n".join(pos_hom_multiple_missense_lowconf),
        "pos_het_multiple_missense_unphaseable": "\n".join(pos_het_multiple_missense_unphaseable),
        "pos_het_multiple_missense_lowconf": "\n".join(pos_het_multiple_missense_lowconf),
        "unknown_var": "\n".join(unknown_var), 
        "n_failed_variants": n_failed,
    }, index=[0])
    report_df.to_csv(args.outfile, index=False, sep='\t')

def main(args):
    generate_variant_report(args)

# EOF
