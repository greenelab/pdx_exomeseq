#!/bin/bash

# Gregory Way 2017
# pdx_exomeseq
# setup-instructions.sh
#
# Setting up environment

# Step 1: Load the python anaconda module
m load python/3.5-Anaconda

# Step 2: Activate the conda environment
conda env create --force --file environment.yml
source activate pdx-exomeseq

# Step 3: Register GATK
# MANUAL STEP: Download Gatk3.8-0 and move to `modules`
# Download from: https://software.broadinstitute.org/gatk/download/
gatk-register modules/GenomeAnalysisTK-3.8-0.tar.bz2

# NOTE: Once this script is run once, to recreate the computational environment, simply run:
# source activate pdx-exomeseq
# However, each time after the conda environment is updated, this script needs to be run again.
