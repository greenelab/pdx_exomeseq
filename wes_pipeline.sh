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
# python scripts/1.run_fastqc.py --command 'fastqc'
#         --data_dir 'data' --output_dir 'results/fastqc_raw' \
#         --walltime '01:30:00' --nodes 1 --cores 4

# Run MultiQC on FastQC results
# multiqc results/fastqc_raw/ \
#         --force --filename 'results/multiqc_report_raw.html'

# Run TrimGalore to cut adapters and filter low quality reads
# TrimGalore is run in `--paired` mode, which performs additional filtering on
# low-quality read pairs between samples. Both pairs of the reads must have
# greater than 20 high quality sequences between them.
# python scripts/2.run_trimgalore.py \
#         --data_dir 'data' --output_dir 'processed/trimmed' \
#         --walltime '04:00:00' --nodes 1 --cores 4

# Run MultiQC again on Trimmed FastQC results from the trimgalore step
# multiqc results/fastqc_trimmed/ \
#         --force --filename 'results/multiqc_report_trimmed.html'

###################
# STEP 2 - Alignment
###################
# Align reads to human reference genome (g1k_v37)
# python scripts/3.run_bwa.py --command 'mem' --genome 'hg' \
#        --data_dir 'processed/trimmed' --output_dir 'processed/sam' \
#        --walltime '05:00:00' --nodes 1 --cores 8

# Also need to align reads to mouse genome (mm9)
# python scripts/3.run_bwa.py --command 'mem' --genome 'mm' \
#        --data_dir 'processed/trimmed' --output_dir 'processed/sam_mouse' \
#        --walltime '05:00:00' --nodes 1 --cores 8

###################
# STEP 3 - Data Conversion and Processing
###################
# Sort SAM and convert to BAM
# python scripts/4.run_samtools.py --command 'sort_name' \
#        --data_dir 'processed/sam' --output_dir 'processed/bam' \
#        --walltime '06:00:00' --nodes 2 --cores 12

# Also need to sort SAM and convert to BAM for mouse
python scripts/4.run_samtools.py --command 'sort_name' \
        --data_dir 'processed/sam_mouse' --output_dir 'processed/bam_mouse' \
        --walltime '06:00:00' --nodes 2 --cores 8

# Here, we need to consider adding a disambiguate step.
# python scripts/6.disambiguate_species.py --command 'disambiguate' \
#        --human_dir 'processed/bam' --mouse_dir 'processed/bam_mouse' \
#        --output_dir 'processed/bam_disambiguate' \
#        --walltime '06:00:00' --nodes 2 --cores 8

# Prep for duplicate removal by cleaning up readpair tags
# python scripts/4.run_samtools.py --command 'fixmate' \
#        --data_dir 'processed/bam' --output_dir 'processed/bam_fixmate' \
#        --walltime '02:30:00' --nodes 2 --cores 4

# Prep for duplicate removal by sorting tagged bam files by position
# python scripts/4.run_samtools.py --command 'sort_position' \
#        --data_dir 'processed/bam_fixmate' \
#        --output_dir 'processed/bam_sort_position' \
#        --walltime '04:30:00' --nodes 2 --cores 8

# Remove duplicate reads
# python scripts/4.run_samtools.py --command 'rmdup' \
#       --data_dir 'processed/bam_sort_position' \
#       --output_dir 'processed/bam_rmdup' \
#       --walltime '02:30:00' --nodes 1 --cores 8

# Create BAM index for duplicate removal
# python scripts/4.run_samtools.py --command 'index_bam' \
#        --data_dir 'processed/bam_rmdup' --output_dir 'processed/bam_rmdup' \
#        --walltime '02:30:00' --nodes 1 --cores 8

###################
# STEP 4 - Variant Calling
###################
# NOTE: Realigning around indels was NOT performed. They are legacy functions
# that are no longer best-practices for GATK HaplotypeCaller pipelines
# (see https://software.broadinstitute.org/gatk/blog?id=7847)

# Assign read groups
# python scripts/5.variant_calling.py --command 'add_read_groups' \
#        --data_dir 'processed/bam_rmdup' --output_dir 'processed/gatk_bam' \
#        --walltime '02:00:00' --nodes 1 --cores 8

# Create index for read group files
# python scripts/5.variant_calling.py --command 'index_bam_gatk' \
#        --data_dir 'processed/gatk_bam' --output_dir 'processed/gatk_bam' \
#        --walltime '2:00:00' --nodes 1 --cores 4

# Call variants with MuTect2
# python scripts/5.variant_calling.py --command 'mutect2' \
#        --data_dir 'processed/gatk_bam' --output_dir 'results/gatk_vcf' \
#        --walltime '05:00:00' --nodes 1 --cores 8

###################
# Step 5 - Filter Mouse Reads from VCF
###################
# python scripts/6.run_mapexr.py --data_dir 'processed/gatk_bam' \
#        --output_dir 'results/mapex_vcf' --vcf_dir 'results/gatk_vcf' \
#        --mapex_path_to_blast_output 'results/blast_out' \
#        --walltime '10:00:00' --nodes 1 --cores 8
