import os

#NGS_PIPELINE_BASE = config['NGS_PIPELINE_BASE']
NGSENV_BASEDIR = os.environ['NGSENV_BASEDIR']
VVG_BASEDIR = os.environ['VVG_BASEDIR']


apptainer_dir = f"{pathlib.Path(VVG_BASEDIR)}/opt/apptainer"
clair3_model_exdir = pathlib.Path(f"{pathlib.Path(VVG_BASEDIR)}/opt/clair3_models").resolve().as_posix()
model_path = config.get('model_path', '/opt/models')
model_name = config.get('model_name', 'r941_prom_sup_g5014')

target_variants_vcf = config.get('target_variants_vcf', '')

rule rename_set_GT:
    input: "{pfx}/{sample}/vcfs/merge_output.vcf.gz"
    output: 
        renamed = temp("{pfx}/{sample}/vcfs/temp.vcf.gz"),
        final = "{pfx}/{sample}/vcfs/variants.vcf.gz",
    params:
        scripts_path = NGSENV_BASEDIR + '/' + 'scripts' + '/' + 'set_gt.py',
        output_dir = lambda w, output: pathlib.Path(output[0]).parent.resolve().as_posix(),
        minimum_depth = config.get('clair3_mindepth', 10),
        minimum_minor_depth = config.get('clair3_minaltdepth', 5),
        minimum_minor_ratio = config.get('clair3_minaltfreq', 0.25),
        headers = f"##cmdline=set_gt.py --minimum_depth {config.get('clair3_mindepth', 10)} --minimum_minor_depth {config.get('clair3_minaltdepth', 5)} --minimum_minor_ratio {config.get('clair3_minaltfreq', 0.25)}"
    shell:
        '''
        bcftools reheader -s {params.output_dir}/sample_id.txt -o {output.renamed} {input} && 
        python3 {params.scripts_path} --infile {output.renamed} --outfile {output.final} --minimum_depth {params.minimum_depth} \
        --minimum_minor_depth {params.minimum_minor_depth} --minimum_minor_ratio {params.minimum_minor_ratio} --headers "{params.headers}"
        '''

rule clair3_out_vcf:
    threads: 4
    input:
        bam = "{pfx}/{sample}/maps/final.bam",
        idx = "{pfx}/{sample}/maps/final.bam.bai"
    output: "{pfx}/{sample}/vcfs/merge_output.vcf.gz"
    params:
        input_dir = lambda w, input: pathlib.Path(input.bam).parent.resolve().as_posix(),
        output_dir = lambda w, output: pathlib.Path(output[0]).parent.resolve().as_posix(),
        ref_dir = pathlib.Path(NGSENV_BASEDIR).as_posix(),
        clair3_extra_flags = config.get('clair3_extra_flags', '') +
            f' --vcf_fn={pathlib.Path(NGSENV_BASEDIR).joinpath(target_variants_vcf).resolve().as_posix()}' if target_variants_vcf else '',
    shell:
        '''
        mkdir -p {params.output_dir} && 
        echo {wildcards.sample} > {params.output_dir}/sample_id.txt && 
        apptainer exec -B {params.input_dir},{params.ref_dir},{params.output_dir},{clair3_model_exdir}:/opt/clair3_models {apptainer_dir}/clair3.sif \
        /opt/bin/run_clair3.sh --bam_fn={params.input_dir}/final.bam --ref_fn={refseq} \
        --threads={threads} --platform=ont --model_path={model_path}/{model_name} --gvcf {params.clair3_extra_flags} \
        --output={params.output_dir}
        '''

