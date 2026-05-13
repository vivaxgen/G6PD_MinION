
# For future reference, when rerio models are available for the latest Clair3 (PyTorch) version

#echo ">>> Cloning ONT rerio repository"
#git clone --depth 1 https://github.com/nanoporetech/rerio.git ${ENVS_DIR}/rerio
#ln -sr ${ENVS_DIR}/rerio ${BASEDIR}/opt/rerio
#ln -sr ${ENVS_DIR}/rerio/clair3_models ${BASEDIR}/opt/clair3_models
#ln -sr ${ENVS_DIR}/G6PD-pipeline/bin/check_n_download_model.py ${BASEDIR}/opt/rerio/check_n_download_model.py

#echo "http://www.bio8.cs.hku.hk/clair3/clair3_models/r941_prom_hac_g360+g422_1235.tar.gz" > ${ENVS_DIR}/rerio/clair3_models/r941_prom_hac_g360+g422_1235_model

#python ${BASEDIR}/opt/rerio/check_n_download_model.py
#mv ${BASEDIR}/opt/clair3_models/ont ${ENVS_DIR}/rerio/clair3_models/r941_prom_hac_g360+g422_1235
#echo "source \${VVG_BASEDIR}/env/G6PD-pipeline/activate.sh" >> ${BASEDIR}/bin/activate.sh
