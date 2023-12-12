
# prepare files for usage

# get the vivaxGEN ngs-pipeline base directory
NGS_PIPELINE_BASE = config['NGS_PIPELINE_BASE']

# include the panel_varcall_lr.smk from vivaxGEN ngs-pipeline
include: f"{NGS_PIPELINE_BASE}/smk/utilities.smk"
include: f"{NGS_PIPELINE_BASE}/smk/general_params.smk"

rule all:
    input:
        f"{refseq}",
        f"{refseq}.fai",
        f"{refseq}.ont.mmi"


# EOF
