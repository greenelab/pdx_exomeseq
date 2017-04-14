#!/bin/bash

# All the steps in order that have been run

# FastQC on all samples
python util/submit_qc_reports.py --data_dir 'data'

# Run TrimGalore to cut adapters and filter low quality reads
python util/run_trimgalore.py

# Run FastQC on the trimmed data
python util/submit_qc_reports.py --data_dir 'data/trimmed'

