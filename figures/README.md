# Main Figures and Legends

Here I outline the main figures the analysis pipeline creates.
I add figure legends as well.

## Read Depth

<object data="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/mosdepth_estimation.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/mosdepth_estimation.pdf">
        This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/mosdepth_estimation.pdf">Download PDF</a>.</p>
    </embed>
</object>

**Figure Legend.** The proportion of the genome covered by thresholded sequencing depth.
Coverage estimation provided by mosdepth ([Pedersen and Quinlan 2018](https://doi.org/10.1093/bioinformatics/btx699)).
The dotted lines represent samples in the lower half of estimated sequencing depth.

## Disambiguating Mouse Reads from PDX models

<object data="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/disambiguate_results.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/disambiguate_results.pdf">
          This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/disambiguate_results.pdf">Download PDF</a>.</p>
    </embed>
</object>

**Figure Legend.** Summary of the proportion of reads assigned to mouse, human, or ambiguously mapped.
Number of read pairs mapped (y axis) plotted against each replicate (x axis) for all samples.
The height of the bars represent the number of reads.
The numbers within each plot represent the proportion of reads mapped to each species (ambiguous; top, human; middle, mouse; bottom).
Disambiguation provided by the disambiguate software ([Ahdesmaki et al. 2017](https://doi.org/10.12688/f1000research.10082.2).

## Filtering and Merging Variants

### Across Replicates

<object data="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/replicates_filtration_results.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/replicates_filtration_results.pdf">
          This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/replicates_filtration_results.pdf">Download PDF</a>.</p>
    </embed>
</object>

**Figure Legend.** Number of variants after each filtration step.
The number of variants filtered is given by the height of the bar _after_ the step is applied.
Each filtration step is applied sequentially.
First, common variation is filtered (gnomAD MAF > 0.05), then low read depth (< 10 reads), then high read depth (> 800 reads).
Each replicate for all samples are shown.

<object data="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/replicates_cosmic_mutcount_results.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/replicates_cosmic_mutcount_results.pdf">
           This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/replicates_cosmic_mutcount_results.pdf">Download PDF</a>.</p>
    </embed>
</object>

**Figure Legend.** Number of COSMIC variants observed in all replicates.
The color of the bars represent the number of _total_ variants in each replicate after filtration.
The text represents the log10 total variant count.

### Merging Replicates

 <object data="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/merged_filtration_results.pdf" type="application/pdf" width="700px" height="700px">
     <embed src="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/merged_filtration_results.pdf">
           This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/merged_filtration_results.pdf">Download PDF</a>.</p>
     </embed>
 </object>

**Figure Legend.** Number of variants after each filtration step after merging replicates.
The number of variants filtered is given by the height of the bar _after_ the step is applied.
Each filtration step is applied sequentially.
First, common variation is filtered (gnomAD MAF > 0.05), then low read depth (<10 reads), then high read depth (> 800 reads).

<object data="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/merged_cosmic_mutcount_results.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/merged_cosmic_mutcount_results.pdf">
           This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/greenelab/pdx_exomeseq/tree/master/figures/merged_cosmic_mutcount_results.pdf">Download PDF</a>.</p>
    </embed>
</object>

**Figure Legend.** Number of COSMIC variants observed after merging replicates.
The color of the bars represent the number of _total_ variants in each replicate after filtration.
The text represents the log10 total variant count.

## Variant Filtration

<object data="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/sift_gnomad/merged_001-F0_sift_gnomad_kde.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/sift_gnomad/merged_001-F0_sift_gnomad_kde.pdf">
           This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/sift_gnomad/merged_001-F0_sift_gnomad_kde.pdf">Download PDF</a>.</p>
    </embed>
</object>

**Figure Legend.** Observations of the effect of filtering variants based on the removal of common variation and read depth.
The intensity of the color represents a higher density of variants.
SIFT score ([Kumar et al. 2009](https://doi.org/10.1038/nprot.2009.86)) represents a prediction of how much the variant is predicted to impact protein function.

**Note:** Plots for all samples (across replicates and merged variants) are located at https://github.com/greenelab/pdx_exomeseq/tree/master/figures/sift_gnomad.

## COSMIC Variants across PDX Passages

<object data="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/upset/upset_sample_001.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/upset/upset_sample_001.pdf">
           This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/upset/upset_sample_001.pdf">Download PDF</a>.</p>
    </embed>
</object>

**Figure Legend.** Exploring the overlap of COSMIC variants across PDX passages.
The UpSet plots ([Conway et al. 2017](https://doi.org/10.1093/bioinformatics/btx364)) depict the number of variants in different sets of the data.
The total overlapping set is highlighted in orange.

A more complicated upset plot is given below:

<object data="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/upset/upset_sample_008.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/upset/upset_sample_008.pdf">
           This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/upset/upset_sample_008.pdf">Download PDF</a>.</p>
    </embed>
</object>

**Figure Legend.** Exploring the overlap of COSMIC variants across PDX passages for sample 008.
The UpSet plots ([Conway et al. 2017](https://doi.org/10.1093/bioinformatics/btx364)) depict the number of variants in different sets of the data.
The total overlapping set is highlighted in orange.

**Note:** Plots for all samples (before and after filtration) are located at https://github.com/greenelab/pdx_exomeseq/tree/master/figures/upset.

## Observing COSMIC Variants across the Population of Samples

### Replicates

<object data="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/cosmic_similarity_replicates.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/cosmic_similarity_replicates.pdf">
           This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/cosmic_similarity_replicates.pdf">Download PDF</a>.</p>
    </embed>
</object>

**Figure Legend.** Displaying the COSMIC profile similarity across replicates.
The Euclidean distance of each COSMIC profile is compared.
High correlation is shown with increasing green intensity.

### After Merging Replicates

<object data="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/cosmic_similarity_merged.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/cosmic_similarity_merged.pdf">
           This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/cosmic_similarity_merged.pdf">Download PDF</a>.</p>
    </embed>
</object>

**Figure Legend.** Displaying the COSMIC profile similarity across samples (replicates are merged).
The Euclidean distance of each COSMIC profile is compared.
High correlation is shown with increasing green intensity.

### Before Variant Filtration (no common variation or depth filtering)

<object data="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/cosmic_prefiltered_similarity_merged.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/cosmic_prefiltered_similarity_merged.pdf">
           This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/cosmic_prefiltered_similarity_merged.pdf">Download PDF</a>.</p>
    </embed>
</object>

**Figure Legend.** Displaying the COSMIC profile similarity across samples before filtration (replicates are merged).
The Euclidean distance of each COSMIC profile is compared.
High correlation is shown with increasing green intensity.

## Oncoprint Diagram of Mutated Genes across Samples

### Replicates

<object data="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/oncoprint_replicates.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/oncoprint_replicates.pdf">
           This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/oncoprint_replicates.pdf">Download PDF</a>.</p>
    </embed>
</object>

**Figure Legend.** Describing mutations in the top 50 most mutated genes by replicate.
A green box indicates the presence of a COSMIC mutation in the given gene for the specific sample.
The top sets of bars represent the number of COSMIC mutations (out of 50) observed in each sample.
The right sets of bars represent the number of samples the gene was observed mutated in (out of 30).
Note that the samples are ordered based on their order in the similarity matrix above.

### Merged Samples

<object data="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/oncoprint_merged.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/oncoprint_merged.pdf">
           This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/greenelab/pdx_exomeseq/blob/master/figures/oncoprint_merged.pdf">Download PDF</a>.</p>
    </embed>
</object>

**Figure Legend.** Describing mutations in the top 50 most mutated genes by sample.
A green box indicates the presence of a COSMIC mutation in the given gene for the specific sample.
The top sets of bars represent the number of COSMIC mutations (out of 50) observed in each sample.
The right sets of bars represent the number of samples the gene was observed mutated in (out of 30).
Note that the samples are ordered based on their order in the similarity matrix above.
