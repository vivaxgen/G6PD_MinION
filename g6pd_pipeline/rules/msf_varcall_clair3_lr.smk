import os

NGS_PIPELINE_BASE = config['NGS_PIPELINE_BASE']
NGSENV_BASEDIR = os.environ['NGSENV_BASEDIR']

apptainer_dir = f"{pathlib.Path(NGS_PIPELINE_BASE).parent.parent}/opt/apptainer"
clair3_model_exdir = pathlib.Path(f"{pathlib.Path(NGS_PIPELINE_BASE).parent.parent}/opt/clair3_models").resolve().as_posix()
model_path = config.get('model_path', '/opt/models')
model_name = config.get('model_name', 'r941_prom_sup_g5014')

target_variants_vcf = config.get('target_variants_vcf', '')

rule clair3_out_vcf:
    threads: 4
    input:
        bam = "{pfx}/{sample}/maps/final.bam",
        idx = "{pfx}/{sample}/maps/final.bam.bai"
    output: "{pfx}/{sample}/vcfs/variants.vcf.gz"
    params:
        input_dir = lambda w, input: pathlib.Path(input.bam).parent.resolve().as_posix(),
        output_dir = lambda w, output: pathlib.Path(output[0]).parent.resolve().as_posix(),
        ref_dir = pathlib.Path(NGSENV_BASEDIR).as_posix(),
        clair3_extra_flags = config.get('clair3_extra_flags', '') +
            f' --vcf_fn={pathlib.Path(NGSENV_BASEDIR).joinpath(target_variants_vcf).resolve().as_posix()}' if target_variants_vcf else '',
    shell:
        'mkdir -p {params.output_dir} && '
        'echo {wildcards.sample} > {params.output_dir}/sample_id.txt && '
        'apptainer exec -B {params.input_dir},{params.ref_dir},{params.output_dir},{clair3_model_exdir}:/opt/clair3_models {apptainer_dir}/clair3.sif'
        ' /opt/bin/run_clair3.sh --bam_fn={params.input_dir}/final.bam --ref_fn={refseq}'
        ' --threads={threads} --platform=ont --model_path={model_path}/{model_name} --gvcf {params.clair3_extra_flags} '
        ' --output={params.output_dir} && '
        'bcftools reheader -s {params.output_dir}/sample_id.txt -o {output} {params.output_dir}/merge_output.vcf.gz'
