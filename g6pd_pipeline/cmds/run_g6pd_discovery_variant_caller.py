
import os
import pathlib
from ngs_pipeline import cerr, cexit, check_NGSENV_BASEDIR
from ngs_pipeline.cmds import run_targeted_variant_caller


snakefiles = {
    'mpileup': 'discovery_mpileup_varcall_lr.smk',
    'freebayes': 'discovery_freebayes_varcall_lr.smk',
    'clair3':  'discovery_clair3_varcall_lr.smk',
}


def init_argparser():
    p = run_targeted_variant_caller.init_argparser()
    p.add_argument('--caller', choices=['mpileup', 'freebayes', 'clair3'],
                   default=None,
                   help='caller to be used [mpileup]')
    p.arg_dict['snakefile'].choices = list(snakefiles.values())
    p.arg_dict['snakefile'].default = snakefiles['mpileup']

    return p


def main(args):

    NGSENV_BASEDIR = pathlib.Path(check_NGSENV_BASEDIR())
    smk_basepath = NGSENV_BASEDIR / 'g6pd_pipeline' / 'rules'

    if args.caller:
        args.snakefile = snakefiles[args.caller]
    args.snakefile = smk_basepath / args.snakefile
    args.target = 'final'
    args.no_config_cascade = True
    args.force = True
    run_targeted_variant_caller.main(args)

# EOF
