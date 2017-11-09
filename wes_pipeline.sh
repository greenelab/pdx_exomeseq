#!/bin/bash

# All the steps in order that have been run

# FastQC on all samples
python scripts/1.run_fastqc.py --data_dir 'data' --output_dir 'results/fastqc_raw' \
        --walltime '01:30:00' --nodes 1 --cores 4

# Run MultiQC on FastQC results
# multiqc results/fastqc_raw

# Run TrimGalore to cut adapters and filter low quality reads
# python scripts/2.run_trimgalore.py --data_dir 'data' --output_dir 'data/trimmed'

# Run FastQC on the trimmed data
# python scripts/1.run_fastqc.py --data_dir 'data/trimmed' --output_dir 'data/fastqc_trimmed'

# Align reads to reference (hg19)
# python scripts/3.run_bwa.py --data_dir 'data/trimmed' --command 'mem'

# Combine paired end reads
# python scripts/3.run_bwa.py --data_dir 'data/trimmed' --command 'sampe'

# Prep to remove duplicate reads
# python scripts/4.run_samtools.py --data_dir 'data/trimmed' --command 'fixmate'

# Prep to remove duplicate reads
# python scripts/4.run_samtools.py --data_dir 'data/trimmed' --command 'sort_order'

# Remove duplicate reads
# python scripts/4.run_samtools.py --data_dir 'data/trimmed' --command 'markdup'

# Adding read groups with picard
# python scripts/5.run_picard.py --data_dir 'data/trimmed' --command 'addreadgroups'
