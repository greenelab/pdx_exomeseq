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

for data_file in wes_files:
    base_name = os.path.basename(data_file)

    command = ['python', os.path.join('util', 'schedule.py'),
               '--command',
               'modules/TrimGalore-0.4.3/trim_galore ' + data_file + ' ' +
               '--path_to_cutadapt /ihome/gway/.conda/envs/pdx-exomeseq/bin/cutadapt ' +
               '--output_dir data/trimmed/ --three_prime_clip_R1 5 --clip_R1 20',
               '--name', 'logs/trimgalore_' + base_name,
               '--walltime', '01:30:00',
               '--filename', 'logs/trimgalore_' + base_name + '.pbs']
    subprocess.call(command)

