#!/bin/bash

# python run_complete_analysis.py --input /storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/regev/matrices/raw/TRNspot_input.h5ad --output regev_test_run --species human --cluster-key merged_response --cluster-key-stratification level_2_annotation --raw-count-layer counts --skip-celloracle --skip-hotspot
python run_complete_analysis.py \
 --input /storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/regev/matrices/raw/TRNspot_input.h5ad \
 --output regev_test_run \
 --species human \
 --cluster-key merged_response \
 --cluster-key-stratification level_2_annotation \
 --raw-count-layer counts \
#  --debug
# python run_complete_analysis.py --input /storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/regev/matrices/raw/TRNspot_input.h5ad --output regev_test_run --species human --cluster-key merged_response --cluster-key-stratification level_2_annotation --raw-count-layer counts
