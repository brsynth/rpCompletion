#!/bin/bash

for channel in `cat ../recipe/conda_channels.txt`; do
  conda config --add channels $channel ;
done
