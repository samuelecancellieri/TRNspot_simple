#!/bin/bash

# RUN TEST ON REGEV DATASET
# python run_complete_analysis.py \
#  --input /storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/regev/matrices/raw/TRNspot_input.h5ad \
#  --output regev_test_run \
#  --species human \
#  --cluster-key merged_response \
#  --cluster-key-stratification level_2_annotation \
#  --raw-count-layer counts

# RUN COMPLETE TEST ON REGEV DATASET
python run_complete_analysis.py \
 --input /storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/regev/matrices/raw/regev_original_with_merged_metadata_scored.h5ad \
 --output regev_complete_run \
 --species human \
 --cluster-key merged_response \
 --clusters "Malignant,Ductal,Acinar,ADM,myCAF,CAF,Macrophage,CD8+ T,CD4+ T,Treg,Natural killer,Adipocyte" \
 --cluster-key-stratification level_2_annotation \
 --raw-count-layer counts
