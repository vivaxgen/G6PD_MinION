
import os
from ngs_pipeline import cerr, cexit
from ngs_pipeline.cmds import run_targeted_variant_caller


def init_argparser():
    return run_targeted_variant_caller.init_argparser()


def main(args):
    args.snakefile = 'panel_varcall_lr.smk'
    args.target = 'merged_report'
    args.no_config_cascade = True
    args.force = True
    run_targeted_variant_caller.main(args)


# EOF

