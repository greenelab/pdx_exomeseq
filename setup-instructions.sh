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

# Step 4: Install custom MAPEXR build from source
# wget --directory-prefix modules/ https://bitbucket.org/bmannakee/mapexr/get/da36687d4585.zip
Rscript util/install_mapexr.R
