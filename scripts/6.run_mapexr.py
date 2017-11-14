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

data_dir = args.data_dir
out_dir = args.output_dir
vcf_dir = args.vcf_dir
blast_out_dir = args.mapex_path_to_blast_output
config = args.config_yaml
walltime = args.walltime
nodes = str(args.nodes)
cores = str(args.cores)

bam_files = []
for path, subdirs, files in os.walk(data_dir):
    for name in files:
        if 'bam.bam_rmdup.bam' in name and '.bam.bai' not in name:
            bam_files.append(name)

# sample_files = [sample_base, sample_bam, sample_bai, sample_vcf]


command_util = os.path.join('util', 'command_wrapper.py')
for sample_1 in bam_files:
    com = ['python', command_util,
           '--sample', sample_1,
           '--command', 'mapex',
           '--mapex_path_to_bam_index', data_dir,
           '--mapex_path_to_vcf', vcf_dir,
           '--mapex_path_to_blast_output', blast_out_dir,
           '--output_directory', out_dir,
           '--config_yaml', config,
           '--walltime', walltime,
           '--nodes', nodes,
           '--cores', cores]
    subprocess.call(com)
