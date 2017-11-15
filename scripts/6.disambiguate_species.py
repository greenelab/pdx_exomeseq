"""
Gregory Way 2017
scripts/6.disambiguate_species.py

Will submit jobs to run variant calling scripts on input bam files
"""

import os
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-h', '--human_dir',
                    help='directory where the human sam files are saved')
parser.add_argument('-m', '--mouse_dir',
                    help='directory where the mouse sam files are saved')
parser.add_argument('-o', '--output_dir',
                    help='Location to save output data',
                    default='.')
parser.add_argument('-y', '--config_yaml',
                    help='Configuration variables for input',
                    default='discovery_variables.yml')
parser.add_argument('-w', '--walltime', default='01:00:00',
                    help='the amount of time alloted to the script')
parser.add_argument('-n', '--nodes', default=1,
                    help='the number of nodes to allocate')
parser.add_argument('-r', '--cores', default=1,
                    help='the number of cores to allocate per node')
args = parser.parse_args()

# Load command arguments
human_dir = args.human_dir
mouse_dir = args.mouse_dir
out_dir = args.output_dir
config = args.config_yaml
walltime = args.walltime
nodes = str(args.nodes)
cores = str(args.cores)

# sample_files will have entry format [sample_id, hg_sam, mm_sam, out_prefix]
sample_files = []
for path, subdirs, files in os.walk(human_dir):
    for name in files:
        sample_id = name.split('.')[0]  # This extracts the sample ID
        human_sam = os.path.join(base_dir, human_dir, name)
        mouse_sam = os.path.join(base_dir, mouse_dir, name)
        output_prefix = os.path.join(base_dir, out_dir, name)
        sample_files.append([sample_id, human_sam, mouse_sam, output_prefix])

sample_files = sample_files[0:2]
command_util = os.path.join('util', 'command_wrapper.py')
for sample_id, human_sam, mouse_sam, output_prefix in sample_files:
    com = ['python', command_util,
           '--sample', sample_id,
           '--command', 'disambiguate',
           '--disambiguate_human', human_sam,
           '--disambiguate_mouse', mouse_sam,
           '--output_directory', output_prefix,
           '--config_yaml', config,
           '--walltime', walltime,
           '--nodes', nodes,
           '--cores', cores]
    print(com)
    subprocess.call(com)
