#!/bin/bash

# All the steps in order that have been run

###################
# STEP 1 - Quality Control
###################
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

###################
# STEP 2 - Alignment
###################
# Align reads to reference (hg19)
# python scripts/3.run_bwa.py --data_dir 'processed/trimmed' --output_dir 'processed/sam' \
#        --command 'mem' --walltime '04:00:00' --nodes 1 --cores 8

###################
# STEP 3 - Data Conversion and Processing
###################
# Sort SAM and convert to BAM
# python scripts/4.run_samtools.py --command 'sort_name' --data_dir 'processed/sam' \
#        --output_dir 'processed/bam' --walltime '06:00:00' --nodes 2 --cores 12

# Prep for duplicate removal by cleaning up readpair tags
# python scripts/4.run_samtools.py --command 'fixmate' --data_dir 'processed/bam' \
#        --output_dir 'processed/bam_fixmate' --walltime '02:30:00' --nodes 2 --cores 4

# Prep for duplicate removal by sorting tagged bam files by position
# python scripts/4.run_samtools.py --command 'sort_position' --data_dir 'processed/bam_fixmate' \
#        --output_dir 'processed/bam_sort_position' --walltime '04:30:00' --nodes 2 --cores 8

# Remove duplicate reads
# python scripts/4.run_samtools.py --command 'rmdup' --data_dir 'processed/bam_sort_position' \
#       --output_dir 'processed/bam_rmdup' --walltime '02:30:00' --nodes 1 --cores 8

# Create BAM index for duplicate removal
# python scripts/4.run_samtools.py --command 'index_bam' --data_dir 'processed/bam_rmdup' \
#        --output_dir 'processed/bam_rmdup' --walltime '02:30:00' --nodes 1 --cores 8

###################
# STEP 4 - Variant Calling
###################
# NOTE: These next two steps were NOT performed. They are legacy functions that are no longer best-practices
# for GATK HaplotypeCaller (see https://software.broadinstitute.org/gatk/blog?id=7847)
# Create fasta index for hg reference
# python util/schedule.py --command 'm load python/3.5-Anaconda && source activate pdx-exomeseq && samtools faidx /lorax/sanchezlab/shared/pdx_exomeseq/reference/human_g1k_v37.fasta' --name 'faindex_hg19' --walltime '03:00:00' --nodes 1 --cores 8 --filename 'logs/faindex_hg19.pbs'

# python scripts/5.variant_calling.py --command 'target_intervals' --data_dir 'processed/bam_rmdup' \
#        --output_dir 'processed/bam_indel_realign' --walltime '03:00:00' --nodes 1 --cores 8

# Haplotype caller but exclude dbSNP vcf
python scripts/5.variant_calling.py --command 'haplotype_caller' --data_dir 'processed/bam_rmdup' \
        --output_dir 'results/gatk_vcf' --walltime '03:00:00' --nodes 1 --cores 8
