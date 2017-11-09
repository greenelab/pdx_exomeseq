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
parser.add_argument('-b', '--sample_2', default=None,
                    help='Second sample of interest')
parser.add_argument('-c', '--command',
                    help='which command to run. Can be one of: \
                         ["fastqc", "trimgalore", "mem", "sampe", "sort", \
                          "markdup"]')
parser.add_argument('-o', '--output_directory',
                    help='the location to save the output files')
parser.add_argument('-y', '--config_yaml', default='discovery_variables.yml',
                    help='Configuration variables for input')
parser.add_argument('-w', '--walltime', default='04:00:00',
                    help='the amount of time alloted to the script')
parser.add_argument('-n', '--nodes', default=1,
                    help='the number of nodes to allocate')
parser.add_argument('-r', '--cores', default=4,
                    help='the number of cores to allocate per node')
args = parser.parse_args()


def schedule_job(command, name, python, nodes=1, cores=4, walltime='04:00:00'):
    output_com = [python, os.path.join('util', 'schedule.py'),
                  '--command', ' '.join(str(x) for x in command),
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
java = config['java']
base_dir = config['directory']
fastqc = config['fastqc']
trimgalore = config['trimgalore']
cutadapt = config['cutadapt']
bwa = config['bwa']
samtools = config['samtools']
picard = config['picard']
gatk = config['gatk']
hg_ref = config['hgreference']

schedule_name = '{}_{}'.format(os.path.basename(sample_1), command)
sample_base = os.path.join(base_dir, output_dir, sample_1.replace('_R1_', '_'))

# Output files
sample_1_sai = sample_base + '_1.sai'
sample_2_sai = sample_base + '_2.sai'
sample_sam = sample_base + '_aln.sam'
sample_sorted_bam = sample_base + '_sorted.bam'
sample_sorted_fixmate_bam = sample_base + '_sorted_fixmate.bam'
sample_sorted_positionsort_bam = sample_base + '_positionsort.bam'
sample_markdup_bam = sample_base + '_markdup.bam'
sample_markdup_bai = sample_base + '_markdup.bam.bai'
sample_readgroup_bam = sample_base + '_markdup.rg.bam'
sample_gatk_intervals = sample_base + '.intervals'
sample_gatk_bam = sample_base + '.GATK.bam'
sample_gatk_bai = sample_base + '.GATK.bam.bai'
sample_gatk_vcf = sample_base + '.GATK.vcf'

# Generate the command calls
fastqc_com = [fastqc, sample_1, '-o', output_dir]
trimgalore_com = [trimgalore, sample_1, '--path_to_cutadapt', cutadapt,
                  '--output_dir', output_dir, 'three_prime_clip_R1', 5,
                  '--clip_R1', 20]
bwa_1_hg_com = [bwa, 'mem', hg_ref, sample_1, '>', sample_1_sai]
bwa_2_hg_com = [bwa, 'mem', hg_ref, sample_2, '>', sample_2_sai]

bwa_paired_com = [bwa, 'sampe', hg_ref, sample_1_sai, sample_2_sai,
                  sample_1, sample_2, '>', sample_sam]

samtools_sort_bam_com = [samtools, 'view', '-@', 4, '-bS', sample_sam, '|',
                         samtools, 'sort', '-n', '-@', 4, '-',
                         '-o', sample_sorted_bam]

samtools_fixmate_com = [samtools, 'fixmate', '-m', sample_sorted_bam,
                        sample_sorted_fixmate_bam]
samtools_positionsort_com = [samtools, 'sort', '-o',
                             sample_sorted_positionsort_bam,
                             sample_sorted_fixmate_bam]
samtools_markdup_com = [samtools, 'markdup', sample_sorted_positionsort_bam,
                        sample_markdup_bam]
samtools_baiindex_com = [samtools, 'index', sample_sorted_positionsort_bam,
                         '-b', sample_markdup_bai]

picard_readgroups_com = [java, '-Xmx50g', '-jar',
                         picard, 'AddOrReplaceReadGroups',
                         'I={}'.format(sample_markdup_bam),
                         'O={}'.format(sample_readgroup_bam),
                         'SORT_ORDER=coordinate',
                         'RGID={}'.format(sample_markdup_bam),
                         'RGLB=bwa-mem', 'RGPL=illumina',
                         'RGSM={}'.format(sample_markdup_bam),
                         'CREATE_INDEX=true', 'RGPU=RGPU',
                         'VALIDATION_STRINGENCY=SILENT']

gatk_localrealign_com = [java, '-Xmx50g', '-jar',
                         gatk, '-T', 'RealignerTargetCreator',
                         '-R', hg_ref, '-I', sample_readgroup_bam,
                         '-O', sample_gatk_intervals]
gatk_realignindels_com = [java, '-Xmx50g', '-jar',
                          gatk, '-allowPotentiallyMisencodedQuals',
                          '-T', 'IndelRealigner', '-R', hg_ref,
                          '-targetIntervals', sample_gatk_intervals,
                          '-I', sample_markdup_bam,
                          '-O', sample_gatk_bam]
gatk_bam_index = [samtools, 'index', sample_gatk_bam, '-b', sample_gatk_bai]
gatk_variant_call = [java, '-Xmx50g', '-jar', gatk, '-T', 'HaplotypeCaller',
                     '-I', sample_gatk_bam, '-o', sample_gatk_vcf,
                     '-R', hg_ref, '--stand_call_conf', 30, '-mmq', 40]

# Schedule a job based on the input command
if command == 'fastqc':
    submit_commands = [fastqc_com]
if command == 'trimgalore':
    submit_commands = [trimgalore_com]
elif command == 'mem':
    submit_commands = [bwa_1_hg_com, bwa_2_hg_com]
elif command == 'sampe':
    submit_commands = [bwa_paired_com]
elif command == 'sort_name':
    submit_commands = [samtools_sort_bam_com]
elif command == 'fixmate':
    submit_commands = [samtools_fixmate_com]
elif command == 'sort_position':
    submit_commands = [samtools_positionsort_com]
elif command == 'markdup':
    submit_commands = [samtools_markdup_com]
elif command == 'index_bam':
    submit_commands = [samtools_baiindex_com]
elif command == 'addreadgroups':
    submit_commands = [picard_readgroups_com]
elif command == 'local_realign':
    submit_commands = [gatk_localrealign_com]
elif command == 'realign_indels':
    submit_commands = [gatk_realignindels_com]
elif command == 'index_gatk':
    submit_commands = [gatk_bam_index]
elif command == 'variant_call_gatk':
    submit_commands = [gatk_variant_call]

if __name__ == '__main__':
    for com in submit_commands:
        schedule_job(command=com, name=schedule_name, python=python,
                     nodes=nodes, cores=cores, walltime=walltime)
