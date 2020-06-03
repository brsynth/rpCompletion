#!/bin/bash

DEBUG='OFF'

size=$1
if [[ "$size" == "" ]]; then
  size="normal"
fi

max_subpaths=$2
if [[ "$max_subpaths" == "" ]]; then
  max_subpaths="10"
fi

python3 $pdb ../src/run.py \
  -rp2paths_compounds in/$size/rp2paths_compounds.csv \
  -rp2_pathways in/$size/rp2_pathways.csv \
  -rp2paths_pathways in/$size/rp2paths_pathways.csv \
  -max_subpaths_filter $max_subpaths \
  -output out/${size}_${max_subpaths} \
  -sm db
