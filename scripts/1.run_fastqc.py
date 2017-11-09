"""
Gregory Way 2017
scripts/1.run_fastqc.py

Call file to perform FastQC on all files each as separate jobs
"""

import os
import argparse
import subprocess
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--data_dir',
                    help='Location to search fmeor files eligible for QC')
parser.add_argument('-o', '--output_dir',
                    help='Location to save output data')
parser.add_argument('-y', '--config_yaml',
                    help='Configuration variables for input',
                    default='discovery_variables.yml')
parser.add_argument('-w', '--walltime', default='04:00:00',
                    help='the amount of time alloted to the script')
parser.add_argument('-n', '--nodes', default=1,
                    help='the number of nodes to allocate')
parser.add_argument('-r', '--cores', default=4,
                    help='the number of cores to allocate per node')
parser.add_argument('-p', '--phenotype', default='pdx_phenotype.csv',
                    help='the filename where the phenotype data is saved')
args = parser.parse_args()

data_dir = args.data_dir
out_dir = args.output_dir
config = args.config_yaml
walltime = args.walltime
nodes = args.nodes
cores = args.cores
phenotype_file = args.phenotype

phenotype_df = pd.read_csv(phenotype_file, index_col=0)

wes_files = []
existing_files = os.listdir(out_dir)
for path, subdirs, files in os.walk(data_dir):
    for name in files:
        if 'fasta.gz' in name or 'fq.gz' in name or 'fastq.gz' in name:
            full_name = os.path.join(path, name)
            fastqc_id = '{}_fastqc.html'.format(name.strip('fastq.gz'))
            if fastqc_id not in existing_files:
                wes_files.append(full_name)

command_util = os.path.join('util', 'command_wrapper.py')
for data_file in wes_files:
    command = ['python', command_util,
               '--sample', data_file,
               '--command', 'fastqc',
               '--output_directory', out_dir,
               '--config_yaml', config,
               '--walltime', walltime,
               '--nodes', nodes,
               '--cores', cores]
    #subprocess.call(command)
