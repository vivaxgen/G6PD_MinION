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
        *[f"{outdir}/{sample}/vcfs/{sample}.clair3.raw.vcf.gz" for sample in IDs],
        f"{outdir}/all.filtered.clair3.vcf.gz",
        #f"{out_dir}/all_f_posf.vcf",
        #f"{out_dir}/all_final.vcf",


#variant calling
apptainer_dir = f"{pathlib.Path(NGS_PIPELINE_BASE).parent.parent}/opt/apptainer"
model_name = 'r941_prom_sup_g5014'

rule clair3:
    threads: 8
    input:
        bam = "{pfx}/{sample}/maps/sorted.bam",
        idx = "{pfx}/{sample}/maps/sorted.bam.bai"
    output:
        vcf = "{pfx}/{sample}/vcfs/{sample}.clair3.raw.vcf.gz"
    log:
        log1 = "{pfx}/{sample}/logs/clair3.log",
    params:
        sample = "{sample}",
        input_dir = lambda w, input: pathlib.Path(input.bam).parent.resolve().as_posix(),
        output_dir = lambda w, output: pathlib.Path(output.vcf).parent.resolve().as_posix(),
        ref_dir = lambda w: pathlib.Path(ngsenv_basedir).parent.as_posix()
    shell:
        'mkdir -p {params.output_dir} && '
        'echo {params.sample} > {params.output_dir}/sample_id.txt && '
        'apptainer exec -B {params.input_dir},{params.ref_dir},{params.output_dir} {apptainer_dir}/clair3.sif'
        ' /opt/bin/run_clair3.sh --bam_fn={params.input_dir}/sorted.bam --ref_fn={refseq}'
        ' --threads={threads} --platform=ont --model_path=/opt/models/{model_name}'
        ' --output={params.output_dir} &&'
        'bcftools reheader -s {params.output_dir}/sample_id.txt -o {output} {params.output_dir}/merge_output.vcf.gz'


#combine raw vcfs

rule combinevcf:
    # aggregation - see the docs here:
    # https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#aggregation
    input:
        vcfs = expand(
            f"{outdir}/{{sample}}/vcfs/{{sample}}.clair3.raw.vcf.gz",
            sample=IDs
        ),
        vcf_idxs = expand(
            f"{outdir}/{{sample}}/vcfs/{{sample}}.clair3.raw.vcf.gz.csi",
            sample=IDs
        )
    output:
        vcf_c = f"{outdir}/all.clair3.raw.vcf.gz"
    log:
        "logs/bcftools-merge.log"
    shell:
        "bcftools merge -o {output.vcf_c} {input.vcfs} 2> {log}"


#filter vcf for quality

rule vcffilter:
    input:
        vcf_c = f"{outdir}/all.clair3.raw.vcf.gz",
    output:
        vcf_f = f"{outdir}/all.filtered.clair3.vcf.gz"
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
