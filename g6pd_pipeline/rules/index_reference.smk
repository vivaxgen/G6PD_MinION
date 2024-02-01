
# prepare files for usage

# get the vivaxGEN ngs-pipeline base directory
NGS_PIPELINE_BASE = config['NGS_PIPELINE_BASE']

# include utilites.smk and general_params.smk from vivaxGEN ngs-pipeline
include: f"{NGS_PIPELINE_BASE}/rules/utilities.smk"
include: f"{NGS_PIPELINE_BASE}/rules/general_params.smk"


rule all:
    input:
        f"{refseq}",
        f"{refseq}.fai",
        f"{refseq}.ont.mmi"


# EOF
