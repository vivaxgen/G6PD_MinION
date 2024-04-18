#!/usr/bin/env bash

$VVGBIN/update-pipeline.sh

python ${VVG_BASEDIR}/opt/rerio/download_model.py --clair3