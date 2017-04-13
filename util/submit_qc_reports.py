"""
Gregory Way 2017
pdx_exomeseq

Call file to perform FastQC on all files each as separate jobs
"""

import os
import subprocess

wes_files = []
for path, subdirs, files in os.walk('data'):
    for name in files:
        wes_files.append(os.path.join(path, name))

for data_file in wes_files:
    base_name = os.path.basename(data_file)
    command = ['python', os.path.join('util', 'schedule.py'),
               '--command', 'modules/FastQC/fastqc ' + data_file,
               '--name', 'fastqc_' + base_name,
               '--walltime', '01:30:00',
               '--filename', 'fastqc_' + base_name + '.pbs']
    subprocess.call(command)

