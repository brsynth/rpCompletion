#!/bin/bash

DEBUG='OFF'

sm=$1
if [[ "$sm" == "" ]]; then
  sm="file"
fi

python3 $pdb ../rpReader/rpReader.py \
  -rp2paths_compounds in/rp2paths_compounds.tsv \
  -rp2_pathways in/rp2_pathways.csv \
  -rp2paths_pathways in/rp2paths_pathways.csv \
  -output out \
  -sm $sm
