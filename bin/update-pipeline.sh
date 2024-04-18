#!/usr/bin/env bash

$VVGBIN/update-pipeline.sh

echo "Checking for new models to download"

python ${VVG_BASEDIR}/opt/rerio/check_n_download_model.py

echo "Update finished"