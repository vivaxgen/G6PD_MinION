__copyright__ = "(C) 2023 Hidayat Trimarsanto, Mariana Barnes"
__license__ = 'MIT'


import os
import pathlib
from glob import glob
from ngs_pipeline import cerr, cexit, check_NGSENV_BASEDIR
from ngs_pipeline.cmds.fetch_clair3_models import fetch_clair3_models, list_available_clair3_models
from ngs_pipeline.cmds import run_targeted_variant_caller
from types import SimpleNamespace

NGSENV_BASEDIR = pathlib.Path(check_NGSENV_BASEDIR())
snakefiles = {
    'freebayes': 'ngs_pipeline::varcaller/freebayes.smk',
    'clair3':  'g6pd_pipeline::msf_varcall_clair3_lr.smk',
}

available_clair3_models = list_available_clair3_models()
available_clair3_models.append("auto_fastq")

def init_argparser():
    p = run_targeted_variant_caller.init_argparser()
    p.add_argument('--caller', choices=['freebayes', 'clair3'],
                   default='freebayes',
                   help='caller to be used [freebayes]')
    p.arg_dict['snakefile'].choices = list(snakefiles.values())
    p.arg_dict['snakefile'].default = snakefiles['freebayes']
    p.add_argument('--clair_model', choices=available_clair3_models,
                   default='auto_fastq',
                   help=('Clair3 models to use, default: [auto_fastq]'
                         'refer to: https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/')
    )
    p.add_argument('--per_amplicon', action='store_true', default=False,
                   help='varcall per amplicon, default: False')
    p.add_argument("--no_flag_failed_variant", action="store_false", default=True,
                   help="Output all variant, including those with depth < mindepth marked with (*) and qual < minqual marked with (^)")
    return p

def get_clair3_path(model):
    model_name_norm = model.strip().replace(".", "").replace("@", "_")
    model_path = f"{os.environ['VVG_BASEDIR']}/opt/clair3-models/{model_name_norm}"
    return model_path, model_name_norm

def ensure_model_exists(model):
    model_path, norm_name = get_clair3_path(model)
    if not os.path.exists(model_path):
        cerr(f"Clair3 model {model} not found at {model_path}. Attempting to fetch...")
        fetch_clair3_models(SimpleNamespace(**{"model": norm_name}))
        if not os.path.exists(model_path):
            cexit(f"Failed to fetch Clair3 model {model}. Please check the model name and your internet connection.")
    return model_path

def check_fastq_for_model(reads):
    from smart_open import open as smart_open
    import re
    known_model = set()
    clair3_models = [m for m in available_clair3_models if m != "auto_fastq"]
    while known_model == set():
        for r in reads:
            with smart_open(r) as f:
                header = f.readline().strip().replace(".", "").replace("@", "_")
                for model in clair3_models:
                    if model in header:
                        known_model.add(model)
    if len(known_model) == 0:
        cexit("Could not determine Clair3 model from FASTQ headers. Please specify the model using --clair_model.")
    return known_model.pop()

def main(args):
    # we will execute targeted variant caller with msf_panel_varcall_lr.smk from vivaxGEN ngs-pipeline
    # see the source here:
    # https://github.com/vivaxgen/ngs-pipeline/blob/main/rules/msf_panel_varcall_lr.smk
    # note: the snakefile is the modular version of panel_varcall_lr.smk
    args.snakefile = 'g6pd_pipeline::g6pd_panel_report.smk'

    # set the target to merged_report
    args.target = ['merge_report', 'all']

    # allow for running outside pipeline base enviroment directory
    args.no_config_cascade = True
    args.force = True
    optional_config = {}

    if args.caller:
        optional_config["varcaller_wf"] = str(snakefiles[args.caller])
        if args.caller == "clair3":
            if args.clair_model == "auto_fastq":
                # Check fastq for model
                resolved_model = check_fastq_for_model(args.infiles)
            else:
                resolved_model = args.clair_model
            optional_config["clair3_model_path"] = ensure_model_exists(resolved_model)
            optional_config["generate_variant_report_extra_flags"] = "--clair3_gvcf"
    if not args.no_flag_failed_variant:
        if "generate_variant_report_extra_flags" in optional_config:
            optional_config["generate_variant_report_extra_flags"] += " --flag_failed_variant"
        optional_config["generate_variant_report_extra_flags"] = "--flag_failed_variant"

    if args.per_amplicon:
        optional_config["amplicon_based"] = True
    run_targeted_variant_caller.run_targeted_variant_caller(args, optional_config)


# EOF
