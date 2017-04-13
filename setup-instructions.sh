#!/bin/bash

# setting up the conda environment is tricky on the cluster

# Step 1: load the python anaconda module
m load python/3.5-Anaconda

# Step 2: activate the conda environment
conda env create --force --file environment.yml
source activate pdx-exomeseq

# Step 3: Unload the anaconda module
m unload python/3.5-Anaconda

# Step 4: Setup channels
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

# Step 5: Install bioconda packages
conda install ngs-disambiguate

