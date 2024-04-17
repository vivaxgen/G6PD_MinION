# Snakefile: discovery_varcall_lr.smk
#
# (C) 2023 Mariana Barnes (mariana.barnes@menzies.edu.au)
#


# get the vivaxGEN ngs-pipeline base directory
NGS_PIPELINE_BASE = config['NGS_PIPELINE_BASE']

# include the panel_varcall_lr.smk from vivaxGEN ngs-pipeline
include: f"{NGS_PIPELINE_BASE}/rules/panel_varcall_lr.smk"

# the following parameters have been defined by panel_varcall_lr.smk above:
# refmap - minimap2 mmi map file
# refseq - fasta file of reference sequence
# target_variants - the file containing list of variants in BED format
# target_regions = the file containing regions in BED format
# variant_info
# infiles - all read files
# outdir
# read_files - a dict with sample codes as keys and read files as values, use this to iterate over samples
# err_files
# min_read_qual
# max_read_len
# min_read_len
# headcrop
# tailcrop


# the following wildcard constraints have been defined by panel_varcall_lr.smk
# as well:
# sample

# the following steps (and their output files) have been defined by rules in
# panel_varcall_lr.smk as well:
# - plot_qc:    {pfx}/reads/raw.fastq.gz >> {pfx}/QC/
# - trim_reads: {pfx}/reads/raw.fastq.gz >> {pfx}/trimmed_reads/trimmed.fastq.gz
# - mapping:    {pfx}/{sample}/trimmed_reads/trimmed.fastq.gz >> {pfx}/{sample}/maps/sorted.bam
#  

# additional variables for this snakemake
min_variant_qual = config.get('min_variant_qual', 30)
min_depth = config.get('min_depth', 10)
IDs = read_files.keys()




rule final:
    input:
        *[f"{outdir}/{sample}/vcfs/{sample}.freebayes.raw.vcf.gz" for sample in IDs],
        f"{outdir}/all.freebayes.filtered.vcf.gz",


# variant calling
rule freebayes_call:
    input:
        bam = "{pfx}/{sample}/maps/sorted.bam",
        idx = "{pfx}/{sample}/maps/sorted.bam.bai"
    output:
        vcf = "{pfx}/{sample}/vcfs/{sample}.freebayes.raw.vcf.gz"
    log:
        log1 = "{pfx}/{sample}/logs/bcftools-freebayes.log",
    shell:
        "freebayes -f {refseq} --target {target_regions} --haplotype-length 0 "
        "--min-base-quality {min_read_qual} {freebayes_extra_flags} {input.bam} 2> log.log1"
        "| bcftools sort -o {output}"


# combine raw vcfs
rule combinevcf:
    # aggregation - see the docs here:
    # https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#aggregation
    input:
        vcfs = expand(
            f"{outdir}/{{sample}}/vcfs/{{sample}}.{{{{caller}}}}.raw.vcf.gz",
            sample=IDs
        ),
        vcf_idxs = expand(
            f"{outdir}/{{sample}}/vcfs/{{sample}}.{{{{caller}}}}.raw.vcf.gz.csi",
            sample=IDs
        )
    output:
        vcf_c = f"{outdir}/all.{{caller}}.raw.vcf.gz"
    log:
        "logs/bcftools-merge-{caller}.log"
    run:
        if len(input.vcfs) == 1:
            shell("ln -P {input.vcfs[0]} {output.vcf_c}")
        else:
            vcfs_to_join = " ".join(input.vcfs)
            shell("bcftools merge -o {output.vcf_c} {vcfs_to_join} 2> {log}")



# filter vcf for quality

rule vcffilter:
    input:
        vcf_c = f"{outdir}/all.{{caller}}.raw.vcf.gz",
    output:
        vcf_f = f"{outdir}/all.{{caller}}.filtered.vcf.gz"
    shell:
        "vcftools --gzvcf {input.vcf_c} --remove-indels --minQ {min_variant_qual} --minDP {min_depth} --recode --stdout | bgzip -c > {output.vcf_f}" 


# run annotation

rule annotate:
    input:
        vcf_p = f"{outdir}/all.{{caller}}.filtered.vcf.gz"
    output:
        vcf_final = f"{outdir}/all.{{caller}}.final.vcf"
    shell:
        "snpEff GRCh38.p14 {input.vcf_p} > {output.vcf_final}"

# EOF
