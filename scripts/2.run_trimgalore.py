"""
Gregory Way 2017
pdx_exomeseq

Call file to perform TrimGalore on all files each as separate jobs
"""

import os
import subprocess

wes_files = []
for path, subdirs, files in os.walk('data'):
    for name in files:
        if 'fastq.gz' in name:
            wes_files.append(os.path.join(path, name))

command_util = os.path.join('util', 'command_wrapper.py')
for data_file in wes_files:
    base_name = os.path.basename(data_file)
    command = ['python', command_util,
               '--sample', data_file,
               '--command', 'trimgalore',
               '--output_dir', os.path.join('data', 'trimmed')]
    subprocess.call(command)
