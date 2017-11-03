"""
Gregory Way 2017
scripts/bwa_mem.py

Will submit jobs to run `bwa mem` on all samples
"""

import os
import subprocess

all_files = []
data_dir = os.path.join('data', 'trimmed')
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
               '--sample_2', sample_2, '--command', 'mem']
    subprocess.call(command)

