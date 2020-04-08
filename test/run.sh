#!/bin/bash

outdir=$1

if [[ "$outdir" == "" ]]; then
  outdir="out/test"
fi

python3 ../src/rpReader.py \
  -rp2paths_compounds in/rp2paths_compounds.tsv \
  -rp2_pathways in/rp2_pathways.csv \
  -rp2paths_pathways in/rp2paths_pathways.csv \
  -output $outdir


  # -upper_flux_bound in/empty_file.csv \
  # -lower_flux_bound 100 \
  # -maxRuleIds 0 \
  # -pathway_id 1000 \
  # -compartment_id 1000 \
  # -species_group_id 1000 \
