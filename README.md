# Whole Exome Sequencing Pipeline for JAX PDX models of Pancreatic Cancer

**Gregory Way<sup>a</sup>, Casey Greene<sup>a</sup>, Yolanda Sanchez<sup>b</sup>**

a. University of Pennsylvania
b. Geisel School of Medicine at Dartmouth

## Summary

Patient derived xenograft (PDX) models were derived from primary and metastatic
tumors from patients admitted to Dartmouth-Hitchcock Medical Center with 
pancreatic adenocarcinoma. The PDX models and tumor samples were whole exome
sequenced (WES) to determine how the mutations from primary tissue and metastases 
propagate or evolve. The following repository outlines the analysis pipeline.

## Pipeline

_This is a work in progress..._

We implemented the [`fastq2vcf`](http://doi.org/10.1186/s13104-015-1027-x) WES processing
pipeline (version 15). The pipeline automatically generages scripts necessary for quality control, 
read filtering, and variant calling. Briefly, `fastq2vcf` first assesses quality using
FastQC, aligns the reads to hg19 using BWA, converts reads to SAM/BAM using samtools,
marks duplicate reads using Picard, realigns using GATK, and then performs variant calling
using four variant callers (GATK UnifiedGenotyper, GATK HaplotypeCaller, SAMtools, and SNVer).

## Compute Environment

All work was performed using the Dartmouth Discovery Cluster Computer

## Steps to Reproduce

```sh
# Install fastq2vcf and dependencies
python util/schedule.py --command "sh install.sh" --name "install_prereq" --walltime "2:59:59" --filename "logs/install.pbs"
```
