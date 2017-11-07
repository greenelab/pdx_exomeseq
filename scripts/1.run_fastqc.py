"""
Gregory Way 2017
pdx_exomeseq

Call file to perform FastQC on all files each as separate jobs
"""

import os
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--data_dir',
                    help='Location to search for files eligible for QC')
parser.add_arguemtn('-o', '--output_dir',
                    help='Location to save fastqc reports')
args = parser.parse_args()

data_dir = args.data_dir
out_dir = args.output_dir

wes_files = []
for path, subdirs, files in os.walk(data_dir):
    for name in files:
        if 'fasta.gz' in name or 'fq.gz' in name:
            wes_files.append(os.path.join(path, name))

command_util = os.path.join('util', 'command_wrapper.py')
for data_file in wes_files:
    base_name = os.path.basename(data_file)
    command = ['python', command_util, '--sample', base_name,
               '--command' 'fastqc',
               '--output_directory', out_dir]
    subprocess.call(command)
