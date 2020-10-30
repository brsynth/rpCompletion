#!/bin/bash

# Pass the engine to be processed by check, if empty all modes will be processed
if [[ $# -gt 0 ]]; then
  mod=$1
else
  mod=$@
fi

bash ./docker/scripts/_indocker.sh flake$mod
bash ./docker/scripts/_indocker.sh bandit
