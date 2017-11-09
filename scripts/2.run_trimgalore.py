"""
Gregory Way 2017
scripts/2.run_trimgalore.py

Call file to perform TrimGalore on all files each as separate jobs
"""

import os
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--data_dir',
                    help='Location of input data')
parser.add_argument('-o', '--output_dir',
                    help='Location to save output data')
parser.add_argument('-y', '--config_yaml',
                    help='Configuration variables for input',
                    default='discovery_variables.yml')
parser.add_argument('-w', '--walltime',
                    help='the amount of time alloted to the script')
parser.add_argument('-n', '--nodes', default=1,
                    help='the number of nodes to allocate')
parser.add_argument('-r', '--cores', default=4,
                    help='the number of cores to allocate per node')
args = parser.parse_args()

data_dir = args.data_dir
out_dir = args.output_dir
config = args.config_yaml
walltime = args.walltime
nodes = args.nodes
cores = args.cores

wes_files = []
for path, subdirs, files in os.walk(data_dir):
    for name in files:
        if 'fastq.gz' in name:
            wes_files.append(os.path.join(path, name))

command_util = os.path.join('util', 'command_wrapper.py')
for data_file in wes_files:
    command = ['python', command_util, '--sample', data_file, '--command',
               'trimgalore', '--output_dir', out_dir,
               '--config_yaml', config, '--walltime', walltime,
               '--nodes', nodes, '--cores', cores]
    subprocess.call(command)
