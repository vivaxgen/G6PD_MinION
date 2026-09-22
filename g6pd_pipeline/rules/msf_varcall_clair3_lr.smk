import os
import ngs_pipeline.rules
include: pkg("ngs_pipeline::varcaller/clair3.smk")

#NGS_PIPELINE_BASE = config['NGS_PIPELINE_BASE']
NGSENV_BASEDIR = os.environ['NGSENV_BASEDIR']

ruleorder: rename_set_GT > clair3_symlink > index_tbi
ruleorder: clair3_g6pd > clair3

use rule clair3 as clair3_g6pd with:
    params:
        contig = "--contig X",


rule rename_set_GT:
    input:
        "<sp>vcfs/clair3/merge_output.vcf.gz"
    output: 
        final = "<sp>vcfs/variants.vcf.gz",
        tbi = "<sp>vcfs/variants.vcf.gz.tbi",
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
        tabix -f {output.final}
        '''


# EOF
