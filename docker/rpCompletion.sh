#!/bin/bash

data_folder=$1
# db_host=$2
# if [[ "$db_host" != "" ]]; then
#   db_host="--link $db_host:$db_host"
# fi


docker-compose run --rm -v $data_folder:$data_folder -w /home/src --entrypoint="" rpcompletion \
  bash

# docker-compose run --rm rpreader $@
