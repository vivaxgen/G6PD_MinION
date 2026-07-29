
echo -e "\e[32m>>>> Linking resource files\e[0m"
${VVGBIN}/link-resource-files.sh ${ENVS_DIR}/G6PD-pipeline/etc/bashrc.d

echo -e "\e[32m>>>> Reloading profiles\e[0m"
reload_vvg_profiles

ENVS_DIR="${ENVS_DIR:-${VVG_BASEDIR}/envs}"
INST_SCRIPTS_DIR="${ENVS_DIR}/G6PD-pipeline/etc/inst-scripts"

if [[ -z ${VVG_MANIFEST_FILE:-} ]]; then
  echo -e "\e[32m>>>> No manifest file provided, installing dependencies with inst-deps.sh\e[0m"
  source ${INST_SCRIPTS_DIR}/inst-deps.sh
fi

echo -e "\e[32m>>>> Indexing reference files\e[0m"
ngs-pl initialize --target panelseq
ngs-pl initialize --target wgs


# EOF
