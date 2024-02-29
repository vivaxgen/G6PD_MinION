
# this is a wrapper to run snakefile index_reference.smk

import os
import pathlib
from ngs_pipeline import cerr, cexit, get_snakefile_path, check_NGSENV_BASEDIR
from ngs_pipeline.cmds import run_snakefile


def init_argparser():
    p = run_snakefile.init_argparser()
    return p


def main(args):

    import g6pd_pipeline

    # NGSENV_BASEDIR is the base directory of the current pipeline (G6PD)
    # NGSENV_BASEDIR = pathlib.Path(check_NGSENV_BASEDIR())
    # smk_basepath = NGSENV_BASEDIR / 'g6pd_pipeline' / 'rules'

    # args.snakefile = smk_basepath / 'index_reference.smk'
    args.snakefile = get_snakefile_path('index_reference.smk',
                                        from_module=g6pd_pipeline)
    args.no_config_cascade = True
    args.force = True

    run_snakefile.main(args)


# EOF
