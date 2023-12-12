
import os
import pathlib
from ngs_pipeline import cerr, cexit, check_NGSENV_BASEDIR
from ngs_pipeline.cmds import run_targeted_variant_caller


def init_argparser():
    p = run_targeted_variant_caller.init_argparser()
    p.arg_dict['snakefile'].choices = ['discovery_mpileup_varcall_lr.smk',
                                       'discovery_clair3_varcall_lr.smk',
                                       'discovery_freebayes_varcall_lr.smk']
    p.arg_dict['snakefile'].default = 'discovery_mpileup_varcall_lr.smk'

    return p


def main(args):

    NGSENV_BASEDIR = pathlib.Path(check_NGSENV_BASEDIR())
    smk_basepath = NGSENV_BASEDIR / 'g6pd_pipeline' / 'smk'

    args.snakefile = smk_basepath / args.snakefile
    args.target = 'final'
    args.no_config_cascade = True
    args.force = True
    run_targeted_variant_caller.main(args)

# EOF
