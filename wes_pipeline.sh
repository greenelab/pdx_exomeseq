#!/bin/bash

# All the steps in order that have been run

# FastQC on all samples
python scripts/1.run_fastqc.py --data_dir 'data' --output_dir 'data/fastqc_raw'

# Run TrimGalore to cut adapters and filter low quality reads
python scripts/2.run_trimgalore.py --data_dir 'data' --output_dir 'data/trimmed'

# Run FastQC on the trimmed data
python scripts/1.run_fastqc.py --data_dir 'data/trimmed' --output_dir 'data/fastqc_trimmed'

