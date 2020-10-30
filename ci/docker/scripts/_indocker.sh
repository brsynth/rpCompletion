#!/bin/bash

source ../extras/.env
source .env

PACKAGE=$PACKAGE \
HOMEDIR=$HOMEDIR \
CONDA_BLD_PATH=$CONDA_BLD_PATH \
docker-compose \
    -f docker/docker-compose.yml \
    --env-file docker/.env \
    build $1

PACKAGE=$PACKAGE \
HOMEDIR=$HOMEDIR \
CONDA_BLD_PATH=$CONDA_BLD_PATH \
docker-compose \
    -f docker/docker-compose.yml \
    --env-file docker/.env \
  run --rm \
  $@
