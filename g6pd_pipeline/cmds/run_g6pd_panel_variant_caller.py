__copyright__ = "(C) 2023 Hidayat Trimarsanto, Mariana Barnes"
__license__ = 'MIT'


import os
from ngs_pipeline import cerr, cexit
from ngs_pipeline.cmds import run_targeted_variant_caller


def init_argparser():
    return run_targeted_variant_caller.init_argparser()


def main(args):
    # we will execute targeted variant caller with msf_panel_varcall_lr.smk from vivaxGEN ngs-pipeline
    # see the source here:
    # https://github.com/vivaxgen/ngs-pipeline/blob/main/rules/msf_panel_varcall_lr.smk
    # note: the snakefile is the modular version of panel_varcall_lr.smk
    args.snakefile = 'msf_panel_varcall_lr.smk'

    # set the target to merged_report
    args.target = 'merged_report'

    # allow for running outside pipeline base enviroment directory
    args.no_config_cascade = True
    args.force = True

    # run targeted variant caller
    run_targeted_variant_caller.main(args)


# EOF
