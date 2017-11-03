"""
Gregory Way 2017
scripts/bwa_wrapper.py

Will perform BWA for human and mouse for paired end reads,
generate sorted BAM files and then disambiguate the reads
"""

import os
import yaml
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--sample',
                    help='First sample of interest')
parser.add_argument('-b', '--sample_2',
                    help='Second sample of interest',
                    default=None)
parser.add_argument('-c', '--command',
                    help='which command to run. Can be one of: \
                         ["fastqc", "mem", "sampe", "sort"]',
                    default='fastqc')
parser.add_argument('-o', '--output_directory',
                    help='the location to save the output files')
parser.add_argument('-y', '--config_yaml',
                    help='Configuration variables for input',
                    default='discovery_variables.yml')
args = parser.parse_args()


def schedule_job(command, name, python, nodes=1, cores=4, walltime='04:00:00'):
    output_com = [python, os.path.join('util', 'schedule.py'),
                  '--command', command,
                  '--name', name,
                  '--walltime', walltime,
                  '--nodes', str(nodes),
                  '--cores', str(cores),
                  '--filename', os.path.join('logs', '{}.pbs'.format(name))]
    return subprocess.call(output_com)

# Load command arguments
sample_1 = args.sample
sample_2 = args.sample_2
command = args.command
output_dir = args.output_directory
config = args.config_yaml

# Load configuration
with open(config, 'r') as stream:
    config = yaml.load(stream)

# Load constants
python = config['python']
base_dir = config['directory']
fasqc = config['fastqc']
bwa = config['bwa']
samtools = config['samtools']
picard = config['picard']
gatk = config['gatk']
hg_ref = config['hg_reference']

schedule_name = '{}_{}'.format(os.path.basename(sample_1), command)
sample_basename = os.path.join(base_dir, output_dir, sample_1.replace('_R1_', '_'))

# Output files
sample_1_sai = sample_basename + '_1.sai'
sample_2_sai = sample_basename + '_2.sai'
sample_sam = sample_basename + '_aln.sam'
sample_sorted_bam = sample_basename + '_sorted.bam'

# Generate the command calls
fastqc_command = [fastqc, sample_1, '-o', output_dir]
bwa_1_hg_com = [bwa, 'mem', hg_ref, sample_1, '>', sample_1_sai]
bwa_2_hg_com = [bwa, 'mem', hg_ref, sample_2, '>', sample_2_sai]

bwa_paired_com = [bwa, 'sampe', hg_ref, sample_1_sai, sample_2_sai,
                  sample_1, sample_2, '>', sample_sam]

samtools_sort_bam_com = [samtools, 'view', '-bS', sample_sam, '|',
                         samtools, 'sort', '-', sample_sorted_bam]

# Schedule a job based on the input command
if command == 'fastqc':
    schedule_job(command=fastqc_command, name=schedule_name, python=python)
elif command == 'mem':
    schedule_job(command=bwa_1_hg_com, name=schedule_name, python=python)
    schedule_job(command=bwa_2_hg_com, name=schedule_name, python=python)
elif command == 'sampe':
    schedule_job(command=bwa_paired_com, name=schedule_name, python=python)
elif command == 'sort':
    schedule_job(command=samtools_sort_bam_com, name=schedule_name, python=python)
