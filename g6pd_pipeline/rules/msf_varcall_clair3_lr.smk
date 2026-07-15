import os
import ngs_pipeline.rules
include: ngs_pipeline.rules.path("ssf_varcall_clair3.smk")

#NGS_PIPELINE_BASE = config['NGS_PIPELINE_BASE']
NGSENV_BASEDIR = os.environ['NGSENV_BASEDIR']

target_variants_vcf = config.get('target_variants_vcf', '')

rule rename_set_GT:
    input: "{pfx}/{sample}/vcfs/merge_output.vcf.gz"
    output: 
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
        python3 {params.scripts_path} --infile {input} --outfile {output.final} --minimum_depth {params.minimum_depth} \
        --minimum_minor_depth {params.minimum_minor_depth} --minimum_minor_ratio {params.minimum_minor_ratio} --headers "{params.headers}"
        '''

use rule ssf_varcall_clair3 as clair3_out_vcf with:
    threads: 4
    input:
        bam = "{pfx}/{sample}/maps/final.bam",
        idx = "{pfx}/{sample}/maps/final.bam.bai",
    output:
        vcf = "{pfx}/{sample}/vcfs/merge_output.vcf.gz",
        vcf_tbi = "{pfx}/{sample}/vcfs/merge_output.vcf.gz.tbi"
    log:
        log1 = "{pfx}/{sample}/logs/clair3.log",
        log2 = "{pfx}/{sample}/logs/clair3.err",
    params:
        sample = lambda wildcards: wildcards.sample,
        platform = config.get('clair3_platform', 'ont'),
        flags = config.get('clair3_flags', '--gvcf'),
        extra_flags = config.get('clair3_extra_flags', '') +
            f' --vcf_fn={pathlib.Path(NGSENV_BASEDIR).joinpath(target_variants_vcf).resolve().as_posix()}' if target_variants_vcf else '',
        model_path = config.get('clair3_model_path', ''), ## fix
        outdir = subpath(output.vcf, parent=True),
        outfmt = ""

