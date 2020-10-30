#!/bin/bash

source .env

bash ./scripts/conda-add_channels.sh

conda build --no-test --output-folder ${CONDA_BLD_PATH} ../recipe
