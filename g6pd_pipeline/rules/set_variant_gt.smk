import os

NGSENV_BASEDIR = os.environ['NGSENV_BASEDIR']

def normalized_depth_file(wildcards):
    if not config.get("neg_control"):
        return []
    return checkpoints.determine_highest_negative_depth.get().output[0]

def normalized_minimum_depth(wildcards, input):
    negative_depth_file = input.get("negative_depth", None)
    if not negative_depth_file:
        return config.get("report_calling_mindepth", 10)
    with open(negative_depth_file) as depth_file:
        return int(depth_file.read().strip())

def normalized_headers(wildcards, input):
    minimum_depth = normalized_minimum_depth(wildcards, input)
    return f"##cmdline=set_gt.py --minimum_depth {minimum_depth}\
        --minimum_minor_depth {config.get('report_calling_minaltdepth', 5)} --minimum_minor_ratio {config.get('report_calling_minaltfreq', 0.25)}\
        {'--drop_zero_depth' if config.get('report_calling_drop_zero_depth', True) else ''}"

rule set_gt:
    input:
        vcf = "<sp>vcfs/variants.vcf.gz",
        idx = "<sp>vcfs/variants.vcf.gz.tbi",
        negative_depth = [],
    output:
        final = "<sp>vcfs/variants.setgt.vcf.gz",
    params:
        scripts_path = NGSENV_BASEDIR + '/' + 'scripts' + '/' + 'set_gt.py',
        output_dir = lambda w, output: pathlib.Path(output[0]).parent.resolve().as_posix(),
        minimum_depth = config.get('report_calling_mindepth', 10),
        minimum_minor_depth = config.get('report_calling_minaltdepth', 5),
        minimum_minor_ratio = config.get('report_calling_minaltfreq', 0.25),
        drop_zero_depth = "--drop_zero_depth" if config.get('report_calling_drop_zero_depth', True) else "",
        headers = f"##cmdline=set_gt.py --minimum_depth {config.get('report_calling_mindepth', 10)}\
        --minimum_minor_depth {config.get('report_calling_minaltdepth', 5)} --minimum_minor_ratio {config.get('report_calling_minaltfreq', 0.25)}\
        {'--drop_zero_depth' if config.get('report_calling_drop_zero_depth', True) else ''}"
    shell:
        '''
        python3 {params.scripts_path} --infile {input.vcf} --outfile {output.final} --minimum_depth {params.minimum_depth} \
        --minimum_minor_depth {params.minimum_minor_depth} --minimum_minor_ratio {params.minimum_minor_ratio} {params.drop_zero_depth} --headers "{params.headers}"
        '''

use rule set_gt as set_gt_norm with:
    input:
        vcf = "<sp>vcfs/variants.vcf.gz",
        idx = "<sp>vcfs/variants.vcf.gz.tbi",
        negative_depth = normalized_depth_file,
    output:
        final = "<sp>vcfs/variants.setgt.norm.vcf.gz",
    params:
        minimum_depth = normalized_minimum_depth,
        headers = normalized_headers,

