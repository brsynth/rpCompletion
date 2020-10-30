#!/bin/bash

source .env
source .secrets
source ../extras/.env

anaconda \
    --token ${ANACONDA_TOKEN} \
    upload \
    --user ${ANACONDA_USER} \
    --label ${ANACONDA_LABEL} \
    ${CONDA_BLD_PATH}/*/${PACKAGE}-*.tar.bz2
