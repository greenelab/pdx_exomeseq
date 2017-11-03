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
args = parser.parse_args()
data_dir = args.data_dir

wes_files = []
for path, subdirs, files in os.walk(data_dir):
    for name in files:
        if 'fasta.gz' in name or 'fq.gz' in name:
            wes_files.append(os.path.join(path, name))

command_util = os.path.join('util', 'command_wrapper.py')
for data_file in wes_files:
    base_name = os.path.basename(data_file)
    command = ['python', command_util, '--sample_1', base_name,
               '--command' 'fastqc']
    subprocess.call(command)

