# Download hg19 genome and generate bwa index files
wget --directory-prefix reference/ ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta.gz
gzip --decompress reference/ucsc.hg19.fasta.gz

# Download 2bit hg19 genome
wget --directory-prefix reference http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

# Download mm9 genome
wget --directory-prefix reference http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit
wget --directory-prefix tools http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa

# Convert 2bit file to fa
cd tools
./twoBitToFa ../reference/hg19.2bit ../reference/hg19.fa
./twoBitToFa ../reference/mm9.2bit ../reference/mm9.fa
cd ..

# Download g1k_v37 (hg19) human fasta file
wget --directory-prefix reference ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
gzip --decompress reference/human_g1k_v37.fasta.gz

# Create combined mouse/human blast database
cp reference/human_g1k_v37.fasta reference/combined_g1k_v37_mm9.fasta
cat reference/mm9.fa | sed 's/chr/mchr/' >> reference/combined_g1k_v37_mm9.fasta
module load blast+/2.6.0
makeblastdb -in reference/combined_g1k_v37_mm9.fasta -parse_seqids -dbtype nucl

# Install Burrows Wheeler Aligner and make hg19 and mm9 index files
python util/schedule.py \
        --command '/ihome/gway/.conda/envs/pdx-exomeseq/bin/bwa index -a bwtsw "reference/human_g1k_v37.fasta"' \
        --name 'hg19-bwa' --walltime '06:00:00' --filename 'bwa-hg19-index.pbs'
python util/schedule.py \
        --command '/ihome/gway/.conda/envs/pdx-exomeseq/bin/bwa index -a bwtsw "reference/mm9.fa"' \
        --name 'mm9-bwa' --walltime '06:00:00' --filename 'bwa-mm9-index.pbs'

# DEPENDENCIES

# GATK
# NOTE: GATK is included in the anaconda environment, but requires a license agreement to actually use.
# MANUAL STEP: Download Gatk3.8-0 https://software.broadinstitute.org/gatk/download/ and move it to `modules` then run:
# gatk-register modules/GenomeAnalysisTK-3.8-0.tar.bz2

wget --directory-prefix reference ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/GATK/All_20170710.vcf.gz
# SAMTOOLS
# Included in the conda environment

# FastQC
# Included in the conda environment

# picard
# included in the conda environment

# vcftools
#wget --directory-prefix modules/ http://pkgs.fedoraproject.org/repo/pkgs/vcftools/vcftools_0.1.11.tar.gz/ddb49e9fa2bfedae36b4dce163adfaa8/vcftools_0.1.11.tar.gz
#tar -xvf modules/vcftools_0.1.11.tar.gz -C modules/ && make --directory modules/vcftools_0.1.11

# tabix
#wget --directory-prefix modules/ https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2
#tar -vxjf modules/tabix-0.2.6.tar.bz2 -C modules/ && make --directory modules/tabix-0.2.6

# SNVer
#wget --directory-prefix modules/ https://sourceforge.net/projects/snver/files/SNVer-0.5.3.tar.gz
#tar -xvf modules/SNVer-0.5.3.tar.gz -C modules/

# MAPEXR
wget --directory-prefix modules/ https://bitbucket.org/bmannakee/mapexr/get/da36687d4585.zip
Rscript util/install_mapexr.R
