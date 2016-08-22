# Here are instructions to install fastq2vcf and dependencies
#wget https://sourceforge.net/projects/fastq2vcf/files/latest/download/fastq2vcf_v15.zip
#unzip fastq2vcf_v15.zip
#rm fastq2vcf_v15.zip
#rm -r __MAC*

# Download hg19 genome and generate bwa index files
#wget --directory-prefix reference/ ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta.gz
#gzip --decompress reference/ucsc.hg19.fasta.gz

# Install Burrows Wheeler Aligner and make hg19 index files
#wget --directory-prefix modules/ https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.5a.tar.bz2 
#tar -vxjf modules/bwa-0.7.5a.tar.bz2 -C modules/ && make --directory modules/bwa-0.7.5a
modules/bwa-0.7.5a/bwa index -a bwtsw "reference/ucsc.hg19.fasta"

