#!/bin/bash

# All the steps in order that have been run

# Fastqc on all samples
python util/submit_qc_reports.py

# Run TrimGalore to cut adapters and filter low quality reads
# First make sure fastqc is in the PATH
export PATH=/lorax/sanchezlab/shared/pdx_exomeseq/modules/FastQC:$PATH
python util/run_trimgalore.py

