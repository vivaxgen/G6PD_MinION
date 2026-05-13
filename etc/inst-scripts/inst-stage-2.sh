
INST_SCRIPTS_DIR="${ENVS_DIR}/G6PD-pipeline/etc/inst-scripts"

echo ">>> Linking resource files"
${VVGBIN}/link-resource-files.sh ${ENVS_DIR}/G6PD-pipeline/etc/bashrc.d

echo ">>>Reloading profiles"
reload_vvg_profiles

echo ">>> Indexing reference sequence"
ngs-pl index-reference

# EOF
