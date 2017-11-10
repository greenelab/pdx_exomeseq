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
                    help='Location of input data')
parser.add_argument('-o', '--output_dir',
                    help='Location to save output data',
                    default='.')
parser.add_argument('-c', '--command',
                    help='One of the bwa commands')
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
command = args.command
config = args.config_yaml
walltime = args.walltime
nodes = str(args.nodes)
cores = str(args.cores)

all_files = []
for path, subdirs, files in os.walk(data_dir):
    for name in files:
        if 'fq.gz' in name:
            all_files.append(name)

paired_reads = []
for name in all_files:
    if '_R1_' in name:
        read_1 = name
        read_2 = name.replace('_R1_', '_R2_')
        read_2 = read_2.replace('val_1', 'val_2')
        paired_reads.append([read_1, read_2])

command_util = os.path.join('util', 'command_wrapper.py')
for sample_1, sample_2 in paired_reads:
    bwa_command = ['python', command_util, '--output_directory', out_dir,
                   '--sample', sample_1, '--sample_2', sample_2,
                   '--command', command, '--config_yaml', config,
                   '--walltime', walltime, '--nodes', nodes,
                   '--cores', cores]
    subprocess.call(bwa_command)
