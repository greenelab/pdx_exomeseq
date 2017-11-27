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

This is a tumor-only analysis, there were no pooled or patient-matched normal
samples available. The following flowchart summarizes the analysis pipeline.

![pdx wes flowchart](figures/pdx_wes_flowchart.png?raw=true)

**A** describes the technical replicates and data-types available across
tumor and mouse passages. **B** outlines our analysis pipeline from quality
control processing raw reads, to alignment and removal of mouse reads, to
calling and annotating variants.

## Pipeline

_This is a work in progress..._

See `wes_pipeline.sh` for our current pipeline.

## Compute Environment

All work was performed using the Dartmouth Discovery Cluster Computer

## Steps to Reproduce

_Work in progress..._
