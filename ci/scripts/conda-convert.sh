#!/bin/bash

source .env
source ../extras/.env

conda convert \
        --platform osx-64 \
        --platform linux-64 \
        --platform win-64 \
        --output-dir ${CONDA_BLD_PATH} \
        ${CONDA_BLD_PATH}/*/${PACKAGE}-*
