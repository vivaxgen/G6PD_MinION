#!/usr/bin/bash

# installation script for vivaxgen g6pd-minion pipeline [https://github.com/vivaxgen/g6pd-minion]

# optional variable:
# - VVG_BASEDIR
# - PIXI_ENVNAME
# - VVG_EXCLUDE
# - VVG_INCLUDE
# - VVG_NGSPL_REPOURL
# - VVG_G6PD_REPOURL
# - VVG_MANIFEST_FILE

__VERSION__="2026.07.15.01"
echo -e "\e[32m>> vivaxGEN G6PD_MinIon pipeline installation script version: ${__VERSION__}\e[0m"


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
if [ -t 0 ] && [ -z "${VVG_VVG_BASEDIR:-}" ]; then
  printf "Pipeline base directory? [./ont-g6pd-pipeline] "
  read VVG_BASEDIR
fi

# default value
VVG_BASEDIR="${VVG_BASEDIR:-./ont-g6pd-pipeline}"

PIXI_ENVNAME='ONT-G6PD'
VVG_EXCLUDE='gatk4'

echo -e "\e[32m>> Installing vivaxGEN G6PD Pipeline pipeline to ${VVG_BASEDIR} with environment name ${PIXI_ENVNAME}\e[0m"

echo -e "\e[32m>> Installing vivaxGEN NGS-Pipeline\e[0m"
source <(curl -L https://raw.githubusercontent.com/vivaxgen/ngs-pipeline/main/install.sh)

echo -e "\e[32m>> Cloning vivaxGEN G6PD Pipeline pipeline\e[0m"
git clone --depth 1  ${VVG_G6PD_REPOURL:-https://github.com/vivaxgen/G6PD_MinION.git} ${ENVS_DIR}/G6PD-pipeline

source ${ENVS_DIR}/G6PD-pipeline/etc/inst-scripts/inst-stage-2.sh

echo "G6PD Pipeline" >> ${ETC_DIR}/installed-repo.txt

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
