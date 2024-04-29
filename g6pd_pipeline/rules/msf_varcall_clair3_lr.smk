import os

NGS_PIPELINE_BASE = config['NGS_PIPELINE_BASE']
NGSENV_BASEDIR = os.environ['NGSENV_BASEDIR']

apptainer_dir = f"{pathlib.Path(NGS_PIPELINE_BASE).parent.parent}/opt/apptainer"
clair3_model_exdir = pathlib.Path(f"{pathlib.Path(NGS_PIPELINE_BASE).parent.parent}/opt/clair3_models").resolve().as_posix()
model_path = config.get('model_path', '/opt/models')
model_name = config.get('model_name', 'r941_prom_sup_g5014')

rule clair3_out_gvcf:
    threads: 4
    input:
        bam = "{pfx}/{sample}/maps/sorted.bam",
        idx = "{pfx}/{sample}/maps/sorted.bam.bai"
    output: "{pfx}/{sample}/vcfs/{sample}.gvcf.gz"
    params:
        input_dir = lambda w, input: pathlib.Path(input.bam).parent.resolve().as_posix(),
        output_dir = lambda w, output: pathlib.Path(output[0]).parent.resolve().as_posix(),
        ref_dir = lambda w: pathlib.Path(ngsenv_basedir).parent.as_posix(),
        clair3_extra_flags = config.get('clair3_extra_flags', '') +
            f' --bed_fn={pathlib.Path(target_variants).resolve().as_posix()}' if target_variants else '',
    shell:
        'mkdir -p {params.output_dir} && '
        'echo {wildcards.sample} > {params.output_dir}/sample_id.txt && '
        'apptainer exec -B {params.input_dir},{params.ref_dir},{params.output_dir},{clair3_model_exdir}:/opt/clair3_models {apptainer_dir}/clair3.sif'
        ' /opt/bin/run_clair3.sh --bam_fn={params.input_dir}/sorted.bam --ref_fn={refseq}'
        ' --threads={threads} --platform=ont --model_path={model_path}/{model_name} --gvcf {params.clair3_extra_flags} '
        ' --output={params.output_dir} &&'
        'bcftools reheader -s {params.output_dir}/sample_id.txt -o {output} {params.output_dir}/merge_output.gvcf.gz'

rule norm_vcf:
    output: "{pfx}/{sample}/vcfs/{sample}.norm.vcf.gz",
    input: "{pfx}/{sample}/vcfs/{sample}.gvcf.gz",
    shell:
        "bcftools convert --gvcf2vcf -f {refseq} -o {output} {input} "


rule filter_vcf:
    input: 
        norm_vcf = "{pfx}/{sample}/vcfs/{sample}.norm.vcf.gz",
        norm_vcf_index = "{pfx}/{sample}/vcfs/{sample}.norm.vcf.gz.csi",
    output: "{pfx}/{sample}/vcfs/variants.vcf.gz",
    shell:
        "bcftools view -R {target_variants} -o {output} {input.norm_vcf}"