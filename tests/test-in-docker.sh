#!/bin/bash

source ../extras/.env

PACKAGE=$PACKAGE \
HOMEDIR=$HOMEDIR \
docker-compose \
    -f pytest/docker/docker-compose.yml \
    --env-file pytest/docker/.env \
    build

# Pass the files to be processed by pytest, if empty all files will be
if [[ $# -eq 1 ]]; then
  file=$1
else
  file=$@
fi


PACKAGE=$PACKAGE \
HOMEDIR=$HOMEDIR \
docker-compose \
    -f pytest/docker/docker-compose.yml \
    --env-file pytest/docker/.env \
  run --rm \
  --entrypoint="" \
  tests sh -c "/opt/conda/lib/python3.7/site-packages/redis_server/bin/redis-server > /dev/null 2>&1 & pytest --verbose $file"
