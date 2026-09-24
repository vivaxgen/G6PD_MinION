import os

NGSENV_BASEDIR = os.environ['NGSENV_BASEDIR']

rule set_gt:
    input:
        vcf = "<sp>vcfs/variants.vcf.gz",
        idx = "<sp>vcfs/variants.vcf.gz.tbi",
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
