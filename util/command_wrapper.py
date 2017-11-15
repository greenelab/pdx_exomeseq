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
                         ["fastqc", "trimgalore", "mem", "sort", "markdup"]')
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

# Output files
sample_sam = sample_base + '.sam'
sample_sorted_bam = sample_base + '_sorted'
sample_sorted_fixmate_bam = sample_base + '_sorted_fixmate.bam'
sample_sorted_positionsort_bam = sample_base + '_positionsort.bam'
sample_markdup_bam = sample_base + '_rmdup.bam'
sample_markdup_bai = sample_base + '.bai'
sample_indel_intervals = sample_base + '.intervals'
sample_indel_realign = sample_base + 'realigned.bam'
sample_addreadgroup = sample_base + '.rg.bam'
sample_bamindex_gatk = sample_base + '.bai'
sample_gatk_vcf = sample_base + '.GATK.vcf'
sample_mapex_out = sample_base + '.tsv'

############################
# Generate the commands
############################
# General purpose module load of pdx-exome seq conda env
conda_build = ['m', 'load', 'python/3.5-Anaconda', '&&',
               'source', 'activate', conda_env, '&&']

mapex_build = ['m', 'load', 'blast+/2.6.0', '&&']

# FastQC
fastqc_com = [fastqc, sample_1, '-o', output_dir]

# TrimGalore
trimgalore_com = [trimgalore, '--paired', sample_1, sample_2,
                  '--output_dir', output_dir,
                  '--fastqc_args', '"--outdir results/fastqc_trimmed/"']

# BWA mem
if command == 'bwa':
    bwa_mem_com = [bwa, 'mem', '-t', '8', hg_ref,
                   os.path.join('processed', 'trimmed', sample_1),
                   os.path.join('processed', 'trimmed', sample_2), '>',
                   sample_sam]

# samtools sort to bam
# `-n` sorts by name, which is required for fixmate
samtools_sort_bam_com = [samtools, 'view', '-bS',
                         os.path.join('processed', 'sam', sample_1),
                         '|', samtools, 'sort', '-n', '-', sample_sorted_bam]

# samtools create fixmate bam
samtools_fixmate_com = [samtools, 'fixmate',
                        os.path.join('processed', 'bam', sample_1),
                        sample_sorted_fixmate_bam]

# samtools sort fixmated bam by position
samtools_positionsort_com = [samtools, 'sort',
                             os.path.join('processed', 'bam_fixmate',
                                          sample_1),
                             sample_sorted_positionsort_bam]

# samtools remove duplicated reads
samtools_rmdup_com = [samtools, 'rmdup',
                      os.path.join('processed', 'bam_sort_position', sample_1),
                      sample_markdup_bam]

# indel realignment - create indel targets
gatk_realigner_com = [gatk, '-T', 'RealignerTargetCreator', '-R', hg_ref,
                      '-I', os.path.join('processed', 'bam_rmdup', sample_1),
                      '-o', sample_indel_intervals]

# samtools create bai indexing in preparation for variant calling
if command == 'index_bam':
    sample_1_file = os.path.join('processed', 'bam_rmdup', sample_1)
    sample_bai_out_file = sample_markdup_bai
elif command == 'index_bam_gatk':
    sample_1_file = os.path.join('processed', 'gatk_bam', sample_1)
    sample_bai_out_file = sample_bamindex_gatk
else:
    sample_1_file = sample_1
    sample_bai_out_file = 'none'

samtools_baiindex_com = [samtools, 'index', sample_1_file, sample_bai_out_file]

# picard add read groups - required for variant calling
picard_addreadgroup_com = [picard, 'AddOrReplaceReadGroups',
                           'I={}'.format(os.path.join('processed', 'bam_rmdup',
                                                      sample_1)),
                           'O={}'.format(sample_addreadgroup),
                           'RGID={}'.format(sample_1),
                           'RGLB=bwa-mem',
                           'RGPL=illumina',
                           'RGSM={}'.format(sample_1),
                           'RGPU={}'.format(sample_1),
                           'CREATE_INDEX=true',
                           'VALIDATION_STRINGENCY=SILENT']

# call variants using GATK MuTect2
gatk_variant_call = [gatk, '-T', 'MuTect2',
                     '--num_cpu_threads_per_data_thread', '8',
                     '--standard_min_confidence_threshold_for_calling', '20',
                     '--min_base_quality_score', '20',
                     '-I:tumor', os.path.join('processed', 'gatk_bam',
                                              sample_1),
                     '-o', sample_gatk_vcf,
                     '-R', hg_ref]

# Remove mouse reads using MAPEX
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

# Schedule a job based on the input command
if command == 'fastqc':
    submit_commands = [fastqc_com]
if command == 'trimgalore':
    # Extending to ensure cutadapt is in path
    conda_build.extend(trimgalore_com)
    submit_commands = [conda_build]
elif command == 'mem':
    conda_build.extend(bwa_mem_com)
    submit_commands = [conda_build]
elif command == 'sort_name':
    conda_build.extend(samtools_sort_bam_com)
    submit_commands = [conda_build]
elif command == 'fixmate':
    conda_build.extend(samtools_fixmate_com)
    submit_commands = [conda_build]
elif command == 'sort_position':
    conda_build.extend(samtools_positionsort_com)
    submit_commands = [conda_build]
elif command == 'rmdup':
    conda_build.extend(samtools_rmdup_com)
    submit_commands = [conda_build]
elif command == 'index_bam' or command == 'index_bam_gatk':
    conda_build.extend(samtools_baiindex_com)
    submit_commands = [conda_build]
elif command == 'target_intervals':
    conda_build.extend(gatk_realigner_com)
    submit_commands = [conda_build]
elif command == 'add_read_groups':
    conda_build.extend(picard_addreadgroup_com)
    submit_commands = [conda_build]
elif command == 'mutect2':
    conda_build.extend(gatk_variant_call)
    submit_commands = [conda_build]
elif command == 'mapex':
    mapex_build.extend(conda_build)
    mapex_build.extend(mapex_remove_mouse_com)
    submit_commands = [mapex_build]

if __name__ == '__main__':
    for com in submit_commands:
        schedule_job(command=com, name=schedule_name, python=python,
                     nodes=nodes, cores=cores, walltime=walltime)
