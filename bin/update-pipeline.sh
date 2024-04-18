#!/usr/bin/env bash

$VVGBIN/update-pipeline.sh

python ${VVG_BASEDIR}/opt/rerio/check_n_download_model.py