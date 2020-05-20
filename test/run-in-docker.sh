#!/bin/bash

cd ../docker
docker-compose run --rm -v $PWD/../test:/home/test -w /home/test --entrypoint="" rpreader \
  sh -c "./run.sh db"
cd -
