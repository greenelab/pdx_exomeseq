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
                         ["fastqc", "trimgalore", "mem", "sort_name", \
                          "sort_position", "fixmate", "rmdup", "index_bam", \
                          "index_bam_gatk", "add_read_groups", \
                          "mutect2", "mapex"]')
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
parser.add_argument('-p', '--mapex_path_to_bam', default='',
                    help='location of the BAM file to use for mapex')
parser.add_argument('-i', '--mapex_path_to_bam_index', default='',
                    help='location of the BAM file index to use for mapex')
parser.add_argument('-v', '--mapex_path_to_vcf', default='',
                    help='location of the VCF file to use for mapex')
parser.add_argument('-l', '--mapex_path_to_blast_output', default='',
                    help='location of the BLAST output file of mapex run')
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
mapex_bam = args.mapex_path_to_bam
mapex_bam_bai = args.mapex_path_to_bam_index
mapex_vcf = args.mapex_path_to_vcf
mapex_blast_out = args.mapex_path_to_blast_output

# Load configuration
with open(config, 'r') as stream:
    config = yaml.load(stream)

# Load constants
python = config['python']
java = config['java']
rscript = config['r']
conda_env = config['condaenv']
hg_ref = config['hgreference']
combined_ref = config['combinedref']
dbsnp = config['dbsnp']
base_dir = config['directory']
fastqc = config['fastqc']
trimgalore = config['trimgalore']
cutadapt = config['cutadapt']
bwa = config['bwa']
samtools = config['samtools']
picard = config['picard']
gatk = config['gatk']
mapex_blast = config['mapexblast']

schedule_name = '{}_{}'.format(os.path.basename(sample_1), command)
sample_name = sample_1.replace('_R1_', '_').replace('_val_1', '')
sample_base = os.path.join(base_dir, output_dir, sample_name)

# Output file suffixes
sample_sam = '{}.sam'.format(sample_base)
sample_sorted_bam = '{}_sorted'.format(sample_base)
sample_sorted_fixmate_bam = '{}_sorted_fixmate.bam'.format(sample_base)
sample_sorted_positionsort_bam = '{}_positionsort.bam'.format(sample_base)
sample_markdup_bam = '{}_rmdup.bam'.format(sample_base)
sample_markdup_bai = '{}.bai'.format(sample_base)
sample_addreadgroup = '{}.rg.bam'.format(sample_base)
sample_bamindex_gatk = '{}.bai'.format(sample_base)
sample_gatk_vcf = '{}.GATK.vcf'.format(sample_base)
sample_mapex_out = '{}.tsv'.format(sample_base)

############################
# Generate the commands
############################
# General purpose module load of pdx-exome seq conda env
conda_build = ['m', 'load', 'python/3.5-Anaconda', '&&',
               'source', 'activate', conda_env, '&&']
mapex_build = ['m', 'load', 'blast+/2.6.0', '&&']

# FastQC
if command == 'fastqc':
    fastqc_com = [fastqc, sample_1, '-o', output_dir]
    conda_build.extend(fastqc_com)

# TrimGalore
if command == 'trimgalore':
    trimgalore_com = [trimgalore, '--paired', sample_1, sample_2,
                      '--output_dir', output_dir,
                      '--fastqc_args', '"--outdir results/fastqc_trimmed/"']
    conda_build.extend(trimgalore_com)

# BWA mem
if command == 'mem':
    bwa_mem_com = [bwa, 'mem', '-t', '8', hg_ref,
                   os.path.join('processed', 'trimmed', sample_1),
                   os.path.join('processed', 'trimmed', sample_2), '>',
                   sample_sam]
    conda_build.extend(bwa_mem_com)

# samtools sort to bam
if command == 'sort_name':
    # `-n` sorts by name, which is required for fixmate
    samtools_sort_bam_com = [samtools, 'view', '-bS',
                             os.path.join('processed', 'sam', sample_1), '|',
                             samtools, 'sort', '-n', '-', sample_sorted_bam]
    conda_build.extend(samtools_sort_bam_com)

# samtools create fixmate bam
if command == 'fixmate':
    samtools_fixmate_com = [samtools, 'fixmate',
                            os.path.join('processed', 'bam', sample_1),
                            sample_sorted_fixmate_bam]
    conda_build.extend(samtools_fixmate_com)

# samtools sort fixmated bam by position
if command == 'sort_position':
    samtools_positionsort_com = [samtools, 'sort',
                                 os.path.join('processed', 'bam_fixmate',
                                              sample_1),
                                 sample_sorted_positionsort_bam]
    conda_build.extend(samtools_positionsort_com)

# samtools remove duplicated reads
if command == 'rmdup':
    samtools_rmdup_com = [samtools, 'rmdup',
                          os.path.join('processed', 'bam_sort_position',
                                       sample_1),
                          sample_markdup_bam]
    conda_build.extend(samtools_rmdup_com)

# samtools create bai indexing in preparation for variant calling
if command == 'index_bam':
    sample_1_file = os.path.join('processed', 'bam_rmdup', sample_1)
    samtools_baiindex_bam_com = [samtools, 'index', sample_1_file,
                                 sample_markdup_bai]
    conda_build.extend(samtools_baiindex_bam_com)

if command == 'index_bam_gatk':
    sample_1_file = os.path.join('processed', 'gatk_bam', sample_1)
    samtools_baiindex_bam_gatk_com = [samtools, 'index', sample_1_file,
                                      sample_bamindex_gatk]
    conda_build.extend(samtools_baiindex_bam_gatk_com)

# picard add read groups - required for variant calling
if command == 'add_read_groups':
    sample_1_file = os.path.join('processed', 'bam_rmdup', sample_1)
    picard_addreadgroup_com = [picard, 'AddOrReplaceReadGroups',
                               'I={}'.format(sample_1_file),
                               'O={}'.format(sample_addreadgroup),
                               'RGID={}'.format(sample_1),
                               'RGLB=bwa-mem',
                               'RGPL=illumina',
                               'RGSM={}'.format(sample_1),
                               'RGPU={}'.format(sample_1),
                               'CREATE_INDEX=true',
                               'VALIDATION_STRINGENCY=SILENT']
    conda_build.extend(picard_addreadgroup_com)

# call variants using GATK MuTect2
if command == 'mutect2':
    gatk_variant_com = [gatk, '-T', 'MuTect2',
                        '--num_cpu_threads_per_data_thread', '8',
                        '--standard_min_confidence_threshold_for_calling', '20',
                        '--min_base_quality_score', '20',
                        '-I:tumor', os.path.join('processed', 'gatk_bam',
                                                 sample_1),
                         '-o', sample_gatk_vcf,
                         '-R', hg_ref]
    conda_build.extend(gatk_variant_com)

# Remove mouse reads using MAPEX
if command == 'mapex':
    mapex_remove_mouse_com = [rscript, '--vanilla', 'util/mapex_wrapper.R',
                              '--path_to_bam', mapex_bam,
                              '--path_to_bam_index', mapex_bam_bai,
                              '--path_to_vcf', mapex_vcf,
                              '--blast_output', mapex_blast_out,
                              '--blast', mapex_blast,
                              '--results_output', sample_mapex_out,
                              '--blast_db', combined_ref,
                              '--num_threads', 8,
                              '--mapq', 1]
    mapex_build.extend(conda_build)
    mapex_build.extend(mapex_remove_mouse_com)

if __name__ == '__main__':
    if command = 'mapex':
        submit_commands = [mapex_build]
    else:
        submit_commands = [conda_build]
    # Submit jobs to cluster
    for com in submit_commands:
        schedule_job(command=com, name=schedule_name, python=python,
                     nodes=nodes, cores=cores, walltime=walltime)
