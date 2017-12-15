#!/bin/bash

# Gregory Way 2017
# pdx_exomeseq
# wes_pipeline.sh
#
# Pipeline for tumor-only WES variant calling. Each script submits a job to the
# DISCOVERY compute cluster at Dartmouth College and will run the specified
# step for all samples

###################
# STEP 1 - Quality Control
###################
# FastQC
python pdx_exomeseq.py fastqc --input_directory 'data' \
         --output_directory 'results/fastqc_raw' \
         --walltime '01:30:00' --nodes 1 --cores 4

# Run MultiQC on FastQC results
python pdx_exomeseq.py multiqc --input_directory 'results/fastqc_raw/' \
         --html_file 'results/multiqc_report_raw.html'

# Run TrimGalore to cut adapters and filter low quality reads
# TrimGalore is run in `--paired` mode, which performs additional filtering on
# low-quality read pairs between samples. Both pairs of the reads must have
# greater than 20 high quality sequences between them.
python pdx_exomeseq.py trimgalore --input_directory 'data' \
         --output_directory 'processed/trimmed' \
         --fastqc_results_dir 'results/fastqc_trimmed/' \
         --walltime '04:00:00' --nodes 1 --cores 4

# Run MultiQC again on Trimmed FastQC results output from the trimgalore step
python pdx_exomeseq.py multiqc --input_directory 'results/fastqc_trimmed/' \
         --html_file 'results/multiqc_report_trimmed.html'

###################
# STEP 2 - Alignment
###################
# Align reads to human reference genome (g1k_v37)
python pdx_exomeseq.py bwa --genome 'hg' --input_directory 'processed/trimmed' \
        --output_directory 'processed/sam' \
        --walltime '05:00:00' --nodes 1 --cores 8

# We also need to align reads to mouse genome (mm9) to disambiguate species
python pdx_exomeseq.py bwa --genome 'mm' --input_directory 'processed/trimmed' \
        --output_directory 'processed/sam_mouse' \
        --walltime '05:00:00' --nodes 1 --cores 8

# NOTE: to save space, delete intermediate files after sorting and compression

###################
# STEP 3 - Remove mouse reads with ngs-disambiguate
###################
# Sort SAM and convert to BAM
python pdx_exomeseq.py samtools --sub_command 'sort_name' --genome 'hg' \
        --input_directory 'processed/sam' \
        --output_directory 'processed/bam' \
        --walltime '06:00:00' --nodes 2 --cores 12

# Also need to sort SAM and convert to BAM for mouse
python pdx_exomeseq.py samtools --sub_command 'sort_name' --genome 'mm' \
        --input_directory 'processed/sam_mouse' \
        --output_directory 'processed/bam_mouse' \
        --walltime '03:00:00' --nodes 2 --cores 12

# Use ngs_disambiguate to identify high confidence human-based reads
python pdx_exomeseq.py disambiguate \
        --human_dir 'processed/bam' --mouse_dir 'processed/bam_mouse' \
        --output_directory 'processed/bam_disambiguate' \
        --walltime '06:00:00' --nodes 2 --cores 8

###################
# STEP 4 - Data Conversion and Processing
###################
# Prep for duplicate removal by cleaning up readpair tags
python pdx_exomeseq.py samtools --sub_command 'fixmate' \
        --input_directory 'processed/bam_disambiguate' \
        --output_directory 'processed/bam_fixmate' \
        --walltime '02:30:00' --nodes 2 --cores 4

# Prep for duplicate removal by sorting tagged bam files by position
python pdx_exomeseq.py samtools --sub_command 'sort_position' \
        --input_directory 'processed/bam_fixmate' \
        --output_directory 'processed/bam_sort_position' \
        --walltime '02:00:00' --nodes 2 --cores 8

# Remove duplicate reads
python pdx_exomeseq.py samtools --sub_command 'rmdup' \
       --input_directory 'processed/bam_sort_position' \
       --output_directory 'processed/bam_rmdup' \
       --walltime '03:30:00' --nodes 1 --cores 4

# Create BAM index for duplicate removal
python pdx_exomeseq.py samtools --sub_command 'index_bam' \
        --input_directory 'processed/bam_rmdup' \
        --output_directory 'processed/bam_rmdup' \
        --walltime '03:30:00' --nodes 1 --cores 4 

###################
# STEP 5 - Variant Calling
###################
# NOTE: Realigning around indels was NOT performed. They are legacy functions
# that are no longer best-practices for GATK HaplotypeCaller pipelines
# (see https://software.broadinstitute.org/gatk/blog?id=7847)

# Assign read groups
python pdx_exomeseq.py variant --sub_command 'add_read_groups' \
        --input_directory 'processed/bam_rmdup' \
        --output_directory 'processed/gatk_bam' \
        --walltime '03:00:00' --nodes 1 --cores 8

# Create index for read group files
python pdx_exomeseq.py samtools --sub_command 'index_bam_gatk' \
        --input_directory 'processed/gatk_bam' \
        --output_directory 'processed/gatk_bam' \
        --walltime '3:30:00' --nodes 1 --cores 4

# Call variants with MuTect2
python pdx_exomeseq.py variant --sub_command 'mutect2' \
        --input_directory 'processed/gatk_bam' \
        --output_directory 'results/gatk_vcf' \
        --walltime '05:00:00' --nodes 1 --cores 8

###################
# STEP 6 - Annotate Variants
###################
# First, download ANNOVAR and associated databases
# Guide: http://annovar.openbioinformatics.org/en/latest/user-guide/startup/
# Db: https://github.com/WGLab/doc-ANNOVAR/blob/master/user-guide/download.md

# The databases we will use are: 
# refGene,cosmic70,gnomad_exome,dbnsfp30a

# First use `convert2annovar` to convert MuTect2 derived VCF files to annovar compatible files
# and then, use `table_annovar` to add annotations as columns to the converted VCF
# python scripts/7.annotate_variants.py
# First use `convert2annovar` to convert MuTect2 derived VCF files to annovar
# compatible files and then use `table_annovar` to add annotations as columns
# to the converted VCF
python scripts/annotate_variants.py

