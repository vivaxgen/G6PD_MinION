__copyright__ = "(C) 2023 Hidayat Trimarsanto, Mariana Barnes"
__license__ = 'MIT'


import os
import pathlib
from glob import glob
from ngs_pipeline import cerr, cexit, check_NGSENV_BASEDIR
from ngs_pipeline.cmds import run_targeted_variant_caller

NGSENV_BASEDIR = pathlib.Path(check_NGSENV_BASEDIR())
snakefiles = {
    'freebayes': 'msf_varcall_freebayes.smk',
    'clair3':  NGSENV_BASEDIR / 'g6pd_pipeline' / 'rules' / 'msf_varcall_clair3_lr.smk',
}

available_clair3_models = [
    a[:-1] for a in glob("*/",
                        root_dir = pathlib.Path(os.environ['NGSENV_BASEDIR']).parent.parent.joinpath("opt/clair3_models")
                        )]
default_models = ['r941_prom_hac_g238', 'r941_prom_sup_g5014', 'r941_prom_hac_g360+g422',
                 'hifi_revio', 'ilmn','hifi_sequel2']
available_clair3_models.extend(default_models)

def init_argparser():
    p = run_targeted_variant_caller.init_argparser()
    p.add_argument('--caller', choices=['freebayes', 'clair3'],
                   default='freebayes',
                   help='caller to be used [freebayes]')
    p.arg_dict['snakefile'].choices = list(snakefiles.values())
    p.arg_dict['snakefile'].default = snakefiles['freebayes']
    p.add_argument('--clair_model', choices=available_clair3_models,
                   default='r941_prom_sup_g5014',
                   help=('Clair3 models to use, default: [r941_prom_sup_g5014]'
                         'refer to: https://github.com/nanoporetech/rerio and'
                         'https://github.com/HKU-BAL/Clair3')
    )
    return p

def get_clair3_path(model):
    if model in default_models:
        # Model already in apptainer
        return '/opt/models'
    else:
        external_model_path = '/opt/clair3_models'
    return external_model_path


def main(args):
    # we will execute targeted variant caller with msf_panel_varcall_lr.smk from vivaxGEN ngs-pipeline
    # see the source here:
    # https://github.com/vivaxgen/ngs-pipeline/blob/main/rules/msf_panel_varcall_lr.smk
    # note: the snakefile is the modular version of panel_varcall_lr.smk
    args.snakefile = 'msf_panel_varcall_lr.smk'

    # set the target to merged_report
    args.target = ['merged_report', 'all'] #  

    # allow for running outside pipeline base enviroment directory
    args.no_config_cascade = True
    args.force = True
    optional_config = {}
    if args.caller:
        print(args.caller)
        optional_config["msf_varcall_wf"] = snakefiles[args.caller]
        if args.caller == "clair3":
            optional_config["model_path"] = get_clair3_path(args.clair_model)
            optional_config["model_name"] = args.clair_model
            optional_config["generate_variant_report_extra_flags"] = "--clair3_gvcf"

    run_targeted_variant_caller.run_targeted_variant_caller(args, optional_config)


# EOF
