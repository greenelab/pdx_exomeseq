"""
Gregory Way 2017
scripts/util/command_wrapper.py

Function that generates commands for each pipeline step
"""

import os
import yaml

import arguments

# Load command arguments
args = arguments.get_args()
print(args)
command = args.which
print(command)
genome = args.genome
input_dir = args.input_directory
output_dir = args.output_directory
config = args.config_yaml
walltime = args.walltime
nodes = str(args.nodes)
cores = str(args.cores)

# Load configuration
with open(config, 'r') as stream:
    config = yaml.load(stream)

# Load constants
python = config['python']
java = config['java']
rscript = config['r']
conda_env = config['condaenv']
if genome == 'hg':
    genome_ref = config['hgreference']
elif genome == 'mm':
    genome_ref = config['mmreference']
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
disambiguate = config['disambiguate']

############################
# Generate the commands
############################
submit_commands = {}

# General purpose module load of pdx-exome seq conda env
conda_build = ['m', 'load', 'python/3.5-Anaconda', '&&',
               'source', 'activate', conda_env, '&&']

all_samples = args.func(args)

# FastQC
if command == 'fastqc':
    for sample_id in all_samples:
        fastqc_com = conda_build + [fastqc, sample_id, '-o', output_dir]
        submit_commands[sample_id] = fastqc_com

# TrimGalore
if command == 'trimgalore':
    for sample_1, sample_2 in all_samples:
        trimgalore_com = [trimgalore, '--paired', sample_1, sample_2,
                          '--output_dir', output_dir,
                          '--fastqc_args',
                          '"--outdir results/fastqc_trimmed/"']
        trimgalore_com = conda_build + trimgalore_com
        submit_commands[sample_id] = trimgalore_com

# BWA mem
if command == 'bwa':
    for sample_1, sample_2 in all_samples:
        sample_name = sample_1.replace('_R1_', '_').replace('_val_1', '')
        sample_base = os.path.join(base_dir, output_dir, sample_name)

        # Output file suffixes
        sample_sam = '{}.sam'.format(sample_base)

        bwa_mem_com = [bwa, 'mem', '-t', '8', genome_ref,
                       os.path.join('processed', 'trimmed', sample_1),
                       os.path.join('processed', 'trimmed', sample_2), '>',
                       sample_sam]
        bwa_mem_com = conda_build + bwa_mem_com
        submit_commands[sample_id] = trimgalore_com

# samtools sort to bam
if command == 'samtools':
    sub_command = args.sub_command
    for sample_id in all_samples:
        sample_name = sample_id.replace('_R1_', '_').replace('_val_1', '')
        sample_base = os.path.join(base_dir, output_dir, sample_name)
        sample_sort_bam = '{}_sorted'.format(sample_base)
        sample_sorted_fixmate_bam = '{}_sorted_fixmate.bam'.format(sample_base)
        sample_sorted_position_bam = '{}_positionsort.bam'.format(sample_base)
        sample_markdup_bam = '{}_rmdup.bam'.format(sample_base)
        sample_markdup_bai = '{}.bai'.format(sample_base)
        sample_bamindex_gatk = '{}.bai'.format(sample_base)

        if sub_command == 'sort_name':
            # `-n` sorts by name, which is required for fixmate
            if genome == 'hg':
                input_samp = os.path.join('processed', 'sam', sample_id)
            elif genome == 'mm':
                input_samp = os.path.join('processed', 'sam_mouse', sample_id)

            samtools_com = [samtools, 'view', '-bS', input_samp, '|',
                            samtools, 'sort', '-n', '-',
                            sample_sort_bam]
        elif sub_command == 'fixmate':
            samtools_com = [samtools, 'fixmate',
                            os.path.join('processed', 'bam_disambiguate',
                                         sample_id),
                            sample_sorted_fixmate_bam]
        elif sub_command == 'sort_position':
            samtools_com = [samtools, 'sort',
                            os.path.join('processed', 'bam_fixmate',
                                         sample_id),
                            sample_sorted_position_bam]

        elif sub_command == 'rmdup':
            samtools_com = [samtools, 'rmdup',
                            os.path.join('processed', 'bam_sort_position',
                                         sample_id),
                            sample_markdup_bam]

        elif sub_command == 'index_bam':
            sample_file = os.path.join('processed', 'bam_rmdup', sample_id)
            samtools_com = [samtools, 'index', sample_file, sample_markdup_bai]

        elif sub_command == 'index_bam_gatk':
            sample_file = os.path.join('processed', 'gatk_bam', sample_id)
            samtools_com = [samtools, 'index', sample_file,
                            sample_bamindex_gatk]

        samtools_com = conda_build + samtools_com
        submit_commands[sample_id] = samtools_com

# Remove mouse reads using ngs_disambiguate
if command == 'disambiguate':
    human_dir = args.human_dir
    mouse_dir = args.mouse_dir

    for sample_id in all_samples:
        disambiguate_com = [disambiguate, '--prefix', sample_id,
                            '--output-dir', output_dir, '--aligner', 'bwa',
                            '--', human_dir, mouse_dir]
        disambiguate_com = conda_build + disambiguate_com
        submit_commands[sample_id] = disambiguate_com


# call variants using GATK MuTect2
if command == 'variant':
    sub_command = args.sub_command
    num_threads = args.threads
    conf = args.min_confidence
    quality = args.min_base_quality_score

    for sample_id in all_samples:
        sample_name = sample_id.replace('_R1_', '_').replace('_val_1', '')
        sample_base = os.path.join(base_dir, output_dir, sample_name)
        sample_gatk_vcf = '{}.GATK.vcf'.format(sample_base)
        sample_addreadgroup = '{}.rg.bam'.format(sample_base)

        if sub_command == 'mutect2':
            tumor_id = os.path.join('processed', 'gatk_bam', sample_id)

            variant_com = [gatk, '-T', 'MuTect2',
                           '--num_cpu_threads_per_data_thread', num_threads,
                           '--standard_min_confidence_threshold_for_calling',
                           conf,
                           '--min_base_quality_score', quality,
                           '-I:tumor', os.path.join('processed', 'gatk_bam',
                                                    sample_id),
                           '-o', sample_gatk_vcf, '-R', genome_ref]
        elif sub_command == 'add_read_groups':
            tumor_id = os.path.join('processed', 'bam_rmdup', sample_id)
            variant_com = [picard, 'AddOrReplaceReadGroups',
                           'I={}'.format(tumor_id),
                           'O={}'.format(sample_addreadgroup),
                           'RGID={}'.format(sample_id),
                           'RGLB=bwa-mem',
                           'RGPL=illumina',
                           'RGSM={}'.format(sample_id),
                           'RGPU={}'.format(sample_id),
                           'CREATE_INDEX=true',
                           'VALIDATION_STRINGENCY=SILENT']
        variant_com = conda_build + variant_com
        submit_commands[sample_id] = variant_com

if __name__ == '__main__':
    # Submit jobs to cluster
    for sample_id, com in submit_commands.items():
        schedule_id = '{}_{}'.format(sample_id, command)
        print(com)
        #arguments.schedule_job(command=com, name=schedule_id, python=python,
        #                       nodes=nodes, cores=cores, walltime=walltime)
