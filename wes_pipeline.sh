#!/bin/bash

# All the steps in order that have been run

# FastQC on all samples
# python scripts/1.run_fastqc.py --data_dir 'data' --output_dir 'results/fastqc_raw' \
#         --walltime '01:30:00' --nodes 1 --cores 4

# Run MultiQC on FastQC results
# multiqc results/fastqc_raw/ --force --filename 'results/multiqc_report_raw.html'

# Run TrimGalore to cut adapters and filter low quality reads
# TrimGalore is run in `--paired` mode, which performs an additional filtering step
# on low-quality read pairs between samples. Both pairs of the reads must have greater
# than 20 high quality sequences between them.
# python scripts/2.run_trimgalore.py --data_dir 'data' --output_dir 'processed/trimmed' \
#         --walltime '04:00:00' --nodes 1 --cores 4

# Run MultiQC again on Trimmed FastQC results from the trimgalore step
# multiqc results/fastqc_trimmed/ --force --filename 'results/multiqc_report_trimmed.html'

# Align reads to reference (hg19)
python scripts/3.run_bwa.py --data_dir 'processed/trimmed' --output_dir 'processed/sam' \
        --command 'mem' --walltime '04:00:00' --nodes 1 --cores 8

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
