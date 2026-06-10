
import os
import pathlib
from glob import glob
from ngs_pipeline import cerr, cexit, check_NGSENV_BASEDIR
from ngs_pipeline.cmds import run_targeted_variant_caller
from ngs_pipeline.cmds.fetch_clair3_models import fetch_clair3_models, list_available_clair3_models
from types import SimpleNamespace
snakefiles = {
    'mpileup': 'discovery_mpileup_varcall_lr.smk',
    'freebayes': 'discovery_freebayes_varcall_lr.smk',
    'clair3':  'discovery_clair3_varcall_lr.smk',
}

available_clair3_models = list_available_clair3_models()
available_clair3_models.append("auto_fastq")

def init_argparser():
    p = run_targeted_variant_caller.init_argparser()
    p.add_argument('--caller', choices=['mpileup', 'freebayes', 'clair3'],
                   default=None,
                   help='caller to be used [mpileup]')
    p.arg_dict['snakefile'].choices = list(snakefiles.values())
    p.arg_dict['snakefile'].default = snakefiles['mpileup']
    p.add_argument('--clair_model', choices=available_clair3_models,
                   default='auto_fastq',
                   help=('Clair3 models to use, default: [auto_fastq]'
                         'refer to: https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/')
    )

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

    NGSENV_BASEDIR = pathlib.Path(check_NGSENV_BASEDIR())
    smk_basepath = NGSENV_BASEDIR / 'g6pd_pipeline' / 'rules'
    optional_config = {}
    if args.caller:
        args.snakefile = snakefiles[args.caller]
        if args.caller == "clair3":
            if args.clair_model == "auto_fastq":
                # Check fastq for model
                resolved_model = check_fastq_for_model(args.infiles)
            else:
                resolved_model = args.clair_model
            optional_config["clair3_model_path"] = ensure_model_exists(resolved_model)

    args.snakefile = smk_basepath / args.snakefile
    args.target = 'final'
    args.no_config_cascade = True
    args.force = True
    run_targeted_variant_caller.run_targeted_variant_caller(args, optional_config)

# EOF
