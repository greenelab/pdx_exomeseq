#!/bin/bash

# Gregory Way 2017
# pdx_exomeseq
# install.sh
#
# Install dependencies for PDX variant calling pipeline

# Download g1k_v37 (hg19) human fasta file
wget --directory-prefix reference ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
gzip --decompress reference/human_g1k_v37.fasta.gz

# Download mm9 genome
wget --directory-prefix reference http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit

# Convert mm9 2bit file to fa
wget --directory-prefix tools http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
cd tools
./twoBitToFa ../reference/mm9.2bit ../reference/mm9.fa
cd ..

# Create BWA index files for hg19 and mm9
python util/schedule.py \
        --command '/ihome/gway/.conda/envs/pdx-exomeseq/bin/bwa index -a bwtsw "reference/human_g1k_v37.fasta"' \
        --name 'hg19-bwa' --walltime '06:00:00' --filename 'bwa-hg19-index.pbs'
python util/schedule.py \
        --command '/ihome/gway/.conda/envs/pdx-exomeseq/bin/bwa index -a bwtsw "reference/mm9.fa"' \
        --name 'mm9-bwa' --walltime '06:00:00' --filename 'bwa-mm9-index.pbs'

# Create samtools fasta index and picard dictionary for hg reference
python util/schedule.py \
        --command 'm load python/3.5-Anaconda && source activate pdx-exomeseq && samtools faidx /lorax/sanchezlab/shared/pdx_exomeseq/reference/human_g1k_v37.fasta' \
        --name 'faindex_hg19' --walltime '03:00:00' --nodes 1 --cores 8 --filename 'logs/faindex_hg19.pbs'

python util/schedule.py \
        --command 'm load python/3.5-Anaconda && source activate pdx-exomeseq && picard CreateSequenceDictionary R=/lorax/sanchezlab/shared/pdx_exomeseq/reference/human_g1k_v37.fasta O=/lorax/sanchezlab/shared/pdx_exomeseq/reference/human_g1k_v37.dict' \
        --name 'fastadict_hg19' --walltime '3:00:00' --nodes 2 --cores 12 --filename 'logs/fastadict_hg19.pbs'

# Download exon data from UCSC table browser
# Selections:
# Assembly = Feb. 2009 (GRCh37/hg19)
# Track = Ensembl Genes
# Output Format = BED
# Output = Coding Exons

# Remove "chr" from the first column of the exon download
# awk '{gsub("chr", "", $0)}1' exon_download.bed > exon_data.bed
