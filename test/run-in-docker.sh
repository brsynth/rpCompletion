#!/bin/bash

size=$1
max_subpaths=$2

cd ../docker
docker-compose run --rm -v $PWD/../test:/home/test -w /home/test --entrypoint="" rpreader \
  sh -c "./run.sh $size $max_subpaths"
cd -
