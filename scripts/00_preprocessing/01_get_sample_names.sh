#!/bin/bash

# SMRT cell 1
# /tgen_labs/jfryer/cores/tgen/r84132_20250721_215156/1_D01

# SMRT cell 2
# /tgen_labs/jfryer/cores/tgen/r84132_20250726_004614/1_C01
cd ../../

cat KBASE_Key.tsv | grep "SMRT3337" | cut -f1,4,7,9,11 > scripts/00_preprocessing/simple_kbase_name.tsv
