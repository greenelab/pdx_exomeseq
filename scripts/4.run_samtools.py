"""
Gregory Way 2017
scripts/4.rum_samtools.py

Will submit jobs to run `samtools` on all samples
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
                    help='One of the samtools commands')
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
nodes = args.nodes
cores = args.cores

sam_files = []
for path, subdirs, files in os.walk(data_dir):
    for name in files:
        if 'fq.gz.sam' in name:
            sam_files.append(name)

command_util = os.path.join('util', 'command_wrapper.py')
for sample_1 in sam_files:
    command = ['python', command_util, '--sample_1', sample_1,
               '--command', command, '--config_yaml', config,
               '--walltime', walltime, '--nodes', nodes, '--cores', cores]
    subprocess.call(command)
