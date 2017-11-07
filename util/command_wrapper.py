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
                         ["fastqc", "trimgalore", "mem", "sampe", "sort"]',
                    default='fastqc')
parser.add_argument('-o', '--output_directory',
                    help='the location to save the output files')
parser.add_argument('-y', '--config_yaml',
                    help='Configuration variables for input',
                    default='discovery_variables.yml')
parser.add_argument('-w', '--walltime',
                    help='the amount of time alloted to the script')
parser.add_argument('-o', '--nodes', default=1,
                    help='the number of nodes to allocate')
parser.add_argument('-r', '--cores', default=1,
                    help='the number of cores to allocate per node')
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
walltime = args.walltime
nodes = args.nodes
cores = args.cores

# Load configuration
with open(config, 'r') as stream:
    config = yaml.load(stream)

# Load constants
python = config['python']
base_dir = config['directory']
fastqc = config['fastqc']
trimgalore = config['trimgalore']
cutadapt = config['cutadapt']
bwa = config['bwa']
samtools = config['samtools']
picard = config['picard']
gatk = config['gatk']
hg_ref = config['hg_reference']

schedule_name = '{}_{}'.format(os.path.basename(sample_1), command)
sample_base = os.path.join(base_dir, output_dir, sample_1.replace('_R1_', '_'))

# Output files
sample_1_sai = sample_base + '_1.sai'
sample_2_sai = sample_base + '_2.sai'
sample_sam = sample_base + '_aln.sam'
sample_sorted_bam = sample_base + '_sorted.bam'

# Generate the command calls
fastqc_com = [fastqc, sample_1, '-o', output_dir]
trimgalore_com = [trimgalore, sample_1, '--path_to_cutadapt', cutadapt,
                  '--output_dir', output_dir, 'three_prime_clip_R1', 5,
                  '--clip_R1', 20]
bwa_1_hg_com = [bwa, 'mem', hg_ref, sample_1, '>', sample_1_sai]
bwa_2_hg_com = [bwa, 'mem', hg_ref, sample_2, '>', sample_2_sai]

bwa_paired_com = [bwa, 'sampe', hg_ref, sample_1_sai, sample_2_sai,
                  sample_1, sample_2, '>', sample_sam]

samtools_sort_bam_com = [samtools, 'view', '-bS', sample_sam, '|',
                         samtools, 'sort', '-', sample_sorted_bam]

# Schedule a job based on the input command
if command == 'fastqc':
    submit_commands = [fastqc_com]
if command == 'trimgalore':
    submit_commands = [trimgalore_com]
elif command == 'mem':
    submit_commands = [bwa_1_hg_com, bwa_2_hg_com]
elif command == 'sampe':
    submit_commands = [bwa_paired_com]
elif command == 'sort':
    submit_commands = [samtools_sort_bam_com]

if __name__ == '__main__':
    for com in submit_commands:
        schedule_job(command=com, name=schedule_name, python=python,
                     nodes=nodes, cores=cores, walltime=walltime)
