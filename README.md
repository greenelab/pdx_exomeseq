# Whole Exome Sequencing Pipeline for JAX FNA-PDX models of Pancreatic Cancer

**Gregory Way<sup>1</sup>, Casey Greene<sup>1</sup>, Yolanda Sanchez<sup>2</sup>**

1. University of Pennsylvania
2. Geisel School of Medicine at Dartmouth

## Summary

Patient derived xenograft (PDX) models were derived from primary and metastatic
tumors from patients admitted to Dartmouth-Hitchcock Medical Center with 
pancreatic adenocarcinoma. The PDX models and tumor samples were whole exome
sequenced (WES) to determine how the mutations from primary tissue and metastases 
propagate or evolve. The following repository outlines the analysis pipeline.

This is a tumor-only analysis; there were no pooled or patient-matched normal
samples available. The following flowchart summarizes the analysis pipeline.

![pdx wes flowchart](figures/pdx_wes_flowchart.png?raw=true)

**Figure 1A** describes the technical replicates and data-types available across
tumor and mouse passages. **Figure 1B** outlines our analysis pipeline from quality
control processing raw reads, to alignment and removal of mouse reads, to
calling and annotating variants.

## Pipeline

See [`wes_pipeline.sh`](https://github.com/greenelab/pdx_exomeseq/blob/master/wes_pipeline.sh)
for our current variant-calling pipeline for tumor-only WES.

## Compute Environment

All work was performed using the Dartmouth Discovery Cluster Computer with the conda
environment specified in
[`environment.yml`](https://github.com/greenelab/pdx_exomeseq/blob/master/environment.yml).

## Steps to Reproduce

There are 3 major steps this repository provides to get from raw sequencing reads to
annotated variants.

#### 1. Setup reproducible computational environment ([`setup_environment.sh`](https://github.com/greenelab/pdx_exomeseq/blob/master/setup_environment.sh), [`install.sh`](https://github.com/greenelab/pdx_exomeseq/blob/master/install.sh))

```bash
# Setup conda environment
bash setup_environment.sh

# NOTE: run `source activate pdx-exomeseq` at the beginning of each session

# Install dependencies and initialize files
# This includes downloading reference genomes and generating several index files
bash install.sh
```

#### 2. Run data processing pipeline ([`wes_pipeline.sh`](https://github.com/greenelab/pdx_exomeseq/blob/master/wes_pipeline.sh))

```bash
# NOTE: the commands in the following script must be run sequentially
# The script will submit several jobs per specified file that can take upwards of
# 12 hours per sample to run _for each command_. This requires the user to specify
# which command is being run by commenting out all others.
bash wes_pipeline.sh
```

Also note that the configuration file `discovery_variables.yml` includes absolute
paths to each tool or resource. It is sufficient to update this file only if paths
to current tools change.

#### 3. Visualize and summarize results

We use Jupyter notebooks and R scripts to visualize and summarize results.

These analyses use a separate conda environment which is specified in
[`analysis_environment.yml`](https://github.com/greenelab/pdx_exomeseq/blob/master/analysis_environment.yml).

Follow these steps to install and begin using this conda environment:

```bash
# Using conda version >= 4.4
conda env create --force --file analysis_environment.yml
conda activate pdx-exomeseq-analysis
```

The following scripts describe our analysis and results:

| Script | Output |
| :----- | :----- |
| [`disambiguate_reads.ipynb`](https://github.com/greenelab/pdx_exomeseq/blob/master/disambiguate_reads.ipynb) | separating mouse reads results |
| [`filter_variants.ipynb`](https://github.com/greenelab/pdx_exomeseq/blob/master/filter_variants.ipynb) | filtration visualization process and filtered VCFs |
| [`variant_allele_frequency.ipynb`](https://github.com/greenelab/pdx_exomeseq/blob/master/variant_allele_frequency.ipynb) | visualize gnomAD by SIFT scores for replicates and filtered merged files |
| [`scripts/viz/variant_overlaps.R`](https://github.com/greenelab/pdx_exomeseq/blob/master/scripts/viz/variant_overlaps.R) | Venn diagrams of variant calling overlaps across replicates and passages |
| [`scripts/viz/oncoprint_similarity.R`](https://github.com/greenelab/pdx_exomeseq/blob/master/scripts/viz/oncoprint_similarity.R) | Oncoprint diagrams and similarity matrices |

