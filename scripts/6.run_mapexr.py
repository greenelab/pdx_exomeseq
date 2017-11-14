"""
Gregory Way 2017
scripts/6.run_mapexr.py

Will submit jobs to run variant calling scripts on input bam files
"""

import os
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--data_dir',
                    help='Location of input data')
parser.add_argument('-i', '--mapex_path_to_bam_index',
                    help='location of the BAM file index to use for mapex')
parser.add_argument('-v', '--mapex_path_to_vcf',
                    help='location of the VCF file to use for mapex')
parser.add_argument('-o', '--output_dir',
                    help='Location to save output data',
                    default='.')
parser.add_argument('-l', '--mapex_path_to_blast_output',
                    help='location of the BLAST output file of mapex run')
parser.add_argument('-v', '--vcf_dir',
                    help='location of the VCF results directory')
parser.add_argument('-y', '--config_yaml',
                    help='Configuration variables for input',
                    default='discovery_variables.yml')
parser.add_argument('-w', '--walltime',
                    help='the amount of time alloted to the script')
parser.add_argument('-n', '--nodes', default=1,
                    help='the number of nodes to allocate')
parser.add_argument('-r', '--cores', default=1,
                    help='the number of cores to allocate per node')
args = parser.parse_args()

# Load command arguments
data_dir = args.data_dir
out_dir = args.output_dir
vcf_dir = args.vcf_dir
blast_out_dir = args.mapex_path_to_blast_output
config = args.config_yaml
walltime = args.walltime
nodes = str(args.nodes)
cores = str(args.cores)

# Set important suffixes
bam_suffix = '.fq.gz.sam_sorted.bam_sorted_fixmate.bam_positionsort.bam' \
             '.bam_rmdup.bam.rg.bam'
bai_suffix = '{}.bai'.format(bam_suffix)
vcf_suffix = '{}.GATK.vcf'.format(bam_suffix)

# sample_files will have entry format [base, bam, bai, vcf]

sample_files = []
for path, subdirs, files in os.walk(data_dir):
    for name in files:
        if 'bam.rg.bam.bai' not in name:
            base = name.split('.')[0]
            samp_bam = os.path.join(data_dir, '{}{}'.format(base, bam_suffix))
            samp_bai = os.path.join(data_dir, '{}{}'.format(base, bai_suffix))
            samp_vcf = os.path.join(vcf_dir, '{}{}'.format(base, vcf_suffix))
            sample_files.append([base, samp_bam, samp_bai, samp_vcf])

command_util = os.path.join('util', 'command_wrapper.py')
for sample_base, sample_bam, sample_bai, sample_vcf in sample_files:
    com = ['python', command_util,
           '--sample', sample_base,
           '--command', 'mapex',
           '--mapex_path_to_bam', sample_bam,
           '--mapex_path_to_bam_index', sample_bai,
           '--mapex_path_to_vcf', sample_vcf,
           '--mapex_path_to_blast_output', blast_out_dir,
           '--output_directory', out_dir,
           '--config_yaml', config,
           '--walltime', walltime,
           '--nodes', nodes,
           '--cores', cores]
    subprocess.call(com)
