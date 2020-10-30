#!/bin/bash

#set -e

for cmd in clean_env build test convert publish
do
  bash ./docker/scripts/conda-${cmd}-indocker.sh
done
