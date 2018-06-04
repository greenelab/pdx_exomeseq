#!/bin/bash

# The following scripts will reproduce the analysis pipeline for the PDX
# whole exome sequencing project. The pipeline will also output a series of
# html notebooks for easy viewing.

# Exit on error
set -o errexit

execute_time=10000000

# Run all files in order
# Notebook 1 - Visualize read depth across the genome
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=$execute_time \
        --execute 1.read-depth-stats.ipynb

# Notebook 2 - Determine the proportion of reads mapping to mouse and human
jupyter nbconvert --to=html \
       --FilesWriter.build_directory=html \
       --ExecutePreprocessor.kernel_name=python3 \
       --ExecutePreprocessor.timeout=$execute_time \
       --execute 2.disambiguate-reads.ipynb

# Notebook 3 - Observe variant SIFT by gnomAD allele frequency
jupyter nbconvert --to=html \
       --FilesWriter.build_directory=html \
       --ExecutePreprocessor.kernel_name=python3 \
       --ExecutePreprocessor.timeout=$execute_time \
       --execute 3.variant-allele-frequency.ipynb

# Notebook 4 - Process the variant calls to output final VCF files
jupyter nbconvert --to=html \
       --FilesWriter.build_directory=html \
       --ExecutePreprocessor.kernel_name=python3 \
       --ExecutePreprocessor.timeout=$execute_time \
       --execute 4.filter-variants.ipynb

# Notebook 5 - Visualize overlaps across patient sets (mutation passage flow)
Rscript --vanilla scripts/nbconverted/5.upset-plots.r

# Notebook 6 - Generate data to visualize common mutations across samples
jupyter nbconvert --to=html \
       --FilesWriter.build_directory=html \
       --ExecutePreprocessor.kernel_name=python3 \
       --ExecutePreprocessor.timeout=$execute_time \
       --execute 6.generate-oncoprint-data.ipynb

# Notebook 7 - Visualize mutations across samples (oncoprint diagrams)
Rscript --vanilla scripts/nbconverted/7.visualize-oncoprint.r

# Convert notebooks to scripts
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb
