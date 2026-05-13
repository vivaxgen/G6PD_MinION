#!/usr/bin/bash

# installation script for G6PD pipeline

# optional variable:
# - VVG_BASEDIR
# - VVG_EXCLUDE
# - VVG_INCLUDE
# - VVG_G6PD_REPOURL

set -eu

# run the base.sh
# Detect the shell from which the script was called
parent=$(ps -o comm $PPID |tail -1)
parent=${parent#-}  # remove the leading dash that login shells have
case "$parent" in
  # shells supported by `micromamba shell init`
  bash|fish|xonsh|zsh)
    shell=$parent
    ;;
  *)
    # use the login shell (basename of $SHELL) as a fallback
    shell=${SHELL##*/}
    ;;
esac

# Parsing arguments
if [ -t 0 ] && [ -z "${VVG_BASEDIR:-}" ]; then
  printf "Pipeline base directory? [./ont-g6pd-pipeline] "
  read VVG_BASEDIR
fi

# default value
VVG_BASEDIR="${VVG_BASEDIR:-./ont-g6pd-pipeline}"

PIXI_ENVNAME='ONT-G6PD'
VVG_EXCLUDE='GATK4'
echo ">> Installing NGS-Pipeline"
source <(curl -L https://raw.githubusercontent.com/vivaxgen/ngs-pipeline/main/install.sh)

#echo Installing apptainer
#micromamba -y install apptainer -c conda-forge -c bioconda
#micromamba -y install squashfuse -c conda-forge

echo ">> Cloning G6PD pipeline"
git clone -depth 1 ${VVG_G6PD_REPOURL:-https://github.com/vivaxgen/G6PD_MinION.git} ${ENVS_DIR}/G6PD-pipeline

echo ">> Executing G6PD pipeline installation stage 2 script"
source ${ENVS_DIR}/G6PD-pipeline/etc/inst-scripts/inst-stage-2.sh

echo "G6PD-pipeline" >> ${ETC_DIR}/inst-envvars

echo
echo "G6PD pipeline has been successfully installed."
echo "To activate the G6PD pipeline environment, either run the activation script"
echo "to get a new shell:"
echo
echo "    $(realpath "${BINDIR}/activate")"
echo
echo "or source the activation script (eg. inside another script):"
echo
echo "    source $(realpath "${BINDIR}/activate")"
echo

# EOF
