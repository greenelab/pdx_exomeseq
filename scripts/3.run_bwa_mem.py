"""
Gregory Way 2017
scripts/3.run_bwa_mem.py

Will submit jobs to run `bwa mem` on all samples
"""

import os
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--data_dir',
                    help='Location to search for files eligible for QC')
parser.add_arguemtn('-o', '--output_dir',
                    help='Location to save fastqc reports')
parser.add_argument('-y', '--config_yaml',
                    help='Configuration variables for input',
                    default='discovery_variables.yml')
parser.add_argument('-w', '--walltime',
                    help='the amount of time alloted to the script')
parser.add_argument('-o', '--nodes', default=1,
                    help='the number of nodes to allocate')
parser.add_argument('-r', '--cores', default=1,
                    help='the number of cores to allocate per node')
args = parser.parse_args()

data_dir = args.data_dir
out_dir = args.output_dir
config = args.config_yaml
walltime = args.walltime
nodes = args.nodes
cores = args.cores

all_files = []
for path, subdirs, files in os.walk(data_dir):
    for name in files:
        if 'trimmed.fq.gz' in name:
            all_files.append(name)

paired_reads = []
for name in all_files:
    if '_R1_' in name:
        read_1 = name
        read_2 = name.replace('_R1_', '_R2_')
        paired_reads.append([os.path.join(data_dir, read_1),
                             os.patj.join(data_dir, read_2)])

command_util = os.path.join('util', 'command_wrapper.py')
for sample_1, sample_2 in paired_reads:
    command = ['python', command_util, '--sample_1', sample_1,
               '--sample_2', sample_2, '--command', 'mem',
               '--config_yaml', config, '--walltime', walltime,
               '--nodes', nodes, '--cores', cores]
    subprocess.call(command)
