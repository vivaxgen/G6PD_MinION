
import os
import pathlib
from ngs_pipeline import cerr, cexit, check_NGSENV_BASEDIR
from ngs_pipeline.cmds import run_snakefile


def init_argparser():
    p = run_snakefile.init_argparser()
    return p


def main(args):

    NGSENV_BASEDIR = pathlib.Path(check_NGSENV_BASEDIR())
    smk_basepath = NGSENV_BASEDIR / 'g6pd_pipeline' / 'smk'

    args.snakefile = smk_basepath / 'index_reference.smk'
    args.no_config_cascade = True
    args.force = True

    run_snakefile.main(args)


# EOF
