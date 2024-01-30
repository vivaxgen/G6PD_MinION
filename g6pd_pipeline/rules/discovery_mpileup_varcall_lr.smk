# Snakefile: discovery_varcall_lr.smk
#
# (C) 2023 Mariana Barnes (mariana.barnes@menzies.edu.au)
#


# get the vivaxGEN ngs-pipeline base directory
NGS_PIPELINE_BASE = config['NGS_PIPELINE_BASE']

# include the panel_varcall_lr.smk from vivaxGEN ngs-pipeline
include: f"{NGS_PIPELINE_BASE}/smk/panel_varcall_lr.smk"

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


# the following wildcard constraints have been defined by panel_varcall_lr.smk as well:
# sample

# the following steps (and their output files) have been defined by rules in panel_varcall_lr.smk as well:
# - plot_qc:    {pfx}/reads/raw.fastq.gz >> {pfx}/QC/
# - trim_reads: {pfx}/reads/raw.fastq.gz >> {pfx}/trimmed_reads/trimmed.fastq.gz
# - mapping:    {pfx}/{sample}/trimmed_reads/trimmed.fastq.gz >> {pfx}/{sample}/maps/sorted.bam
#  

# additional variables for this snakemake
min_variant_qual = config.get('min_variant_qual', 30)
min_depth = config.get('min_depth', 10)
IDs = read_files.keys()

#def final_filelist(w):
#    return [f"{out_dir}/{sample_id}/{sample_id}.out.txt" for sample_id in IDs]


rule final:
    input:
        *[f"{outdir}/{sample}/vcfs/{sample}.mpileup.raw.vcf.gz" for sample in IDs],
        f"{outdir}/all.filtered.mpileup.vcf.gz",
        #f"{out_dir}/all_f_posf.vcf",
        #f"{out_dir}/all_final.vcf",


# XXX: please consider rule index_ref and index_minimap into a separate, stand-alone shell scripts
# that needs to be run just once during the set up.

# index reference if not already done

#rule index_ref:
#    input:
#        ref_genome
#    output:
#        f"{ref_genome}.fai"
#    shell:
#        "samtools faidx {input}"

#rule index_minimap:
#    input:
#        ref_genome
#    output:
#        f"{ref_genome}.mmi"
#    shell:
#        "minimap2 -x map-ont -d {output} {input}"


#variant calling
rule mpileup_call:
    input:
        bam = "{pfx}/{sample}/maps/sorted.bam",
        idx = "{pfx}/{sample}/maps/sorted.bam.bai"
    output:
        vcf = "{pfx}/{sample}/vcfs/{sample}.mpileup.raw.vcf.gz"
    log:
        log1 = "{pfx}/{sample}/logs/bcftools-mpileup.log",
        log2 = "{pfx}/{sample}/logs/bcftools-call.log"
    shell:
        "bcftools mpileup -f {refseq} -R {target_regions} -a AD,ADF,ADR,DP,SP,SCR {input.bam} 2> {log.log1} "
        "| bcftools call -mv -o {output} 2> {log.log2}"


#combine raw vcfs

rule combinevcf:
    # aggregation - see the docs here:
    # https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#aggregation
    input:
        vcfs = expand(
            f"{outdir}/{{sample}}/vcfs/{{sample}}.mpileup.raw.vcf.gz",
            sample=IDs
        ),
        vcf_idxs = expand(
            f"{outdir}/{{sample}}/vcfs/{{sample}}.mpileup.raw.vcf.gz.csi",
            sample=IDs
        )
    output:
        vcf_c = f"{outdir}/all.mpileup.raw.vcf.gz"
    log:
        "logs/bcftools-merge.log"
    shell:
        "bcftools merge -o {output.vcf_c} {input.vcfs} 2> {log}"


#filter vcf for quality

rule vcffilter:
    input:
        vcf_c = f"{outdir}/all.mpileup.raw.vcf.gz",
    output:
        vcf_f = f"{outdir}/all.filtered.mpileup.vcf.gz"
    shell:
        "vcftools --gzvcf {input.vcf_c} --remove-indels --minQ {min_variant_qual} --minDP {min_depth} --recode --stdout | bgzip -c > {output.vcf_f}" 

#normalizing hets

#rule normhets:
#    input:
#        vcf_f=f"{out_dir}/all.filtered.vcf"
#    output:
#        vcf_fh=f"{out_dir}/all.het_filtered.vcf"
#    shell:


#run annotation

rule annotate:
    input:
        vcf_p = f"{outdir}/all_f_posf.vcf"
    output:
        vcf_final = f"{outdir}/all_final.vcf"
    shell:
        "snpEff GRCh38.p14 {input.vcf_p} > {output.vcf_final}"


# utilities

# EOF
