"""
Gregory Way 2017
scripts/bwa_wrapper.py

Will perform BWA for human and mouse for paired end reads,
generate sorted BAM files and then disambiguate the reads
"""

import os
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--sample_1',
                    help='Sample for paired end read 1')
parser.add_argument('-b', '--sample_2',
                    help='Sample for paired end read 2')
args = parser.parse_args()

sample_1 = args.sample_1
sample_2 = args.sample_2

# Load constants
base_dir = '/lorax/sanchezlab/shared/pdx_exomeseq/'
bwa = base_dir + 'modules/bwa-0.7.5a/bwa'
hg_ref_file = base_dir + 'reference/ucsc.hg19.fasta'
mm_ref_file = base_dir + 'reference/mm9.fa'
sample_basename = base_dir + 'data/bwa/' + sample_1.replace('_R1_', '_')

# Output files
sample_1_hg_sai = sample_basename + '_hg_1.sai'
sample_2_hg_sai = sample_basename + '_hg_2.sai'
sample_1_mm_sai = sample_basename + '_mm_1.sai'
sample_2_mm_sai = sample_basename + '_mm_2.sai'

sample_hg_sam = sample_basename + '_hg_aln.sam'
sample_mm_sam = sample_basename + '_mm_aln.sam'

sample_hg_sorted_bam = sample_basename + '_hg_sorted.bam'
sample_mm_sorted_bam = sample_basename + '_mm_sorted.bam'

disambig_dir = 'data/disambiguated/'

# Generate the command calls
bwa_1_hg_com = [bwa, 'aln', hg_ref_file, sample_1, '>', sample_1_hg_sai]
bwa_2_hg_com = [bwa, 'aln', hg_ref_file, sample_2, '>', sample_2_hg_sai]
bwa_1_mm_com = [bwa, 'aln', mm_ref_file, sample_1, '>', sample_1_mm_sai]
bwa_2_mm_com = [bwa, 'aln', mm_ref_file, sample_2, '>', sample_2_mm_sai]

bwa_paired_hg_com = [bwa, 'sampe', hg_ref_file, sample_1_hg_sai, sample_2_hg_sai,
                     sample_1, sample_2, '>', sample_hg_sam]
bwa_paired_mm_com = [bwa, 'sampe', mm_ref_file, sample_1_mm_sai, sample_2_mm_sai,
                     sample_1, sample_2, '>', sample_mm_sam]

samtools_sorted_bam_hg_com = ['samtools', 'view', '-bS', sample_hg_sam, '|',
                              'samtools', 'sort', '-', sample_hg_sorted_bam]
samtools_sorted_bam_mm_com = ['samtools', 'view', '-bS', sample_mm_sam, '|',
                              'samtools', 'sort', '-', sample_mm_sorted_bam]

disambiguate_com = ['ngs-disambiguate', sample_hg_sorted_bam, sample_mm_sorted_bam,
                    '--aligner', 'bwa', '--no-sort', '--output-dir', disambig_dir]

# Make the commands in sequence
#subprocess.call(bwa_1_hg_com)
subprocess.call(bwa_2_hg_com)
#subprocess.call(bwa_1_mm_com)
#subprocess.call(bwa_2_mm_com)
#subprocess.call(bwa_paired_hg_com)
#subprocess.call(bwa_paired_mm_com)
#subprocess.call(samtools_sorted_bam_hg_com)
#subprocess.call(samtools_sorted_bam_mm_com)
#subprocess.call(disambiguate_com)

