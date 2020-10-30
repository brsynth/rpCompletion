#!/bin/bash

source .env
source ../extras/.env

bash ./scripts/conda-add_channels.sh

conda build --test ${CONDA_BLD_PATH}/*/${PACKAGE}-*.tar.bz2
