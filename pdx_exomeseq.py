"""
Gregory Way 2017
pdx_exomeseq

Main function that aggregates all samples and determines which command to run
depending on command line arguments. Then submits these commands as scheduled
jobs to the Dartmouth Discovery cluster.

Usage:

Invoke via command line with specific arguments including:

{'fastqc', 'multiqc', 'trimgalore', 'bwa', 'samtools', 'disambiguate',
'variant'}

Some commands may also have subcommands. E.g. with samtools, specify one of the
following subcommands:

{'sort_name', 'fixmate', 'sort_position', 'rmdup', 'index_bam',
'index_bam_gatk'}

E.g.

    python pdx_exomeseq.py samtools
            --sub_command 'sort_name' \
            --genome 'hg' \
            --input_directory 'processed/sam' \
            --output_directory 'processed/bam' \
            --walltime '06:00:00' \
            --nodes 2 \
            --cores 12

The specific pipeline using this script is given in `wes_pipeline.sh`
"""

import os
import yaml

import util.arguments as arguments

# Load command arguments
args = arguments.get_args()
command = args.which

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
exonbed = config['exonbed']
base_dir = config['directory']
fastqc = config['fastqc']
multiqc = config['multiqc']
trimgalore = config['trimgalore']
cutadapt = config['cutadapt']
bwa = config['bwa']
samtools = config['samtools']
picard = config['picard']
gatk = config['gatk']
disambiguate = config['disambiguate']
mosdepth = config['mosdepth']

############################
# Generate the commands
############################
submit_commands = {}

# General purpose module load of pdx-exome seq conda env
conda_build = ['m', 'load', 'python/3-Anaconda', '&&',
               'source', 'activate', conda_env, '&&']
java_load = ['m', 'load', 'java/1.8', '&&']

# Obtain samples to process
all_samples = args.func(args)

# FastQC
if command == 'fastqc':
    for sample_id in all_samples:
        fastqc_com = conda_build + [fastqc, sample_id, '-o', output_dir]
        submit_commands[sample_id] = fastqc_com

# MultiQC
if command == 'multiqc':
    html_file = args.html_file

    multiqc_com = conda_build + [multiqc, input_dir, '--force', '--filename',
                                 html_file]
    submit_commands['multiqc'] = multiqc_com

# TrimGalore
if command == 'trimgalore':
    fastqc_results_dir = args.fastqc_results_dir

    for sample_1, sample_2 in all_samples:
        # The Second output specifies which directory to perform fastqc on
        trimgalore_com = [trimgalore,
                          '--paired', sample_1, sample_2,
                          '--output_dir', output_dir,
                          '--fastqc_args',
                          '"--outdir {}"'.format(fastqc_results_dir)]
        trimgalore_com = conda_build + trimgalore_com
        submit_commands[sample_1] = trimgalore_com

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
        submit_commands[sample_sam] = bwa_mem_com

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
            sample_file = os.path.join(input_dir, sample_id)
            sample_output_file = os.path.join(output_dir, '{}.bai'.format(sample_id))
            samtools_com = [samtools, 'index', sample_file, sample_output_file]

        elif sub_command == 'flagstat':
            sample_file = os.path.join(input_dir, sample_id)
            sample_output_file = os.path.join(output_dir, '{}.flagstat.txt'.format(sample_id))
            samtools_com = [samtools, 'flagstat', sample_file, '>', sample_output_file]

        elif sub_command == 'merge':
            sample_file = os.path.join(input_dir, sample_id)
            if 'L001' in sample_id and '.bai' not in sample_id:
                output_bam = '{}.merged.bam'.format(sample_id.split('_')[0])
                output_bam = os.path.join(output_dir, output_bam)

                replicate_2 = sample_file.replace('L001', 'L002')
                replicate_3 = sample_file.replace('L001', 'L003')
                replicate_4 = sample_file.replace('L001', 'L004')
                samtools_com = [samtools, 'merge', '-r', output_bam, sample_file,
                                replicate_2, replicate_3, replicate_4]
            else:
                continue

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
        sample_addreadgroup = '{}.rg.bam'.format(sample_base)

        if sub_command == 'mutect2':
            tumor_id = os.path.join(input_dir, sample_id)
            output_file = os.path.join(output_dir, '{}.GATK.vcf'.format(sample_id))

            variant_com = [gatk, '-T', 'MuTect2',
                           '--num_cpu_threads_per_data_thread', num_threads,
                           '--standard_min_confidence_threshold_for_calling',
                           conf,
                           '--min_base_quality_score', quality,
                           '-I:tumor', tumor_id, '-o', output_file,
                           '-R', genome_ref]

        elif sub_command == 'add_read_groups':
            tumor_id = os.path.join(input_dir, sample_id)
            output_file = os.path.join(output_dir, '{}.rg.bam'.format(sample_id))

            variant_com = [picard, 'AddOrReplaceReadGroups',
                           'I={}'.format(tumor_id),
                           'O={}'.format(output_file),
                           'RGID={}'.format(sample_id),
                           'RGLB=bwa-mem',
                           'RGPL=illumina',
                           'RGSM={}'.format(sample_id),
                           'RGPU={}'.format(sample_id),
                           'CREATE_INDEX=true',
                           'VALIDATION_STRINGENCY=SILENT']

        variant_com = conda_build + java_load + variant_com
        submit_commands[sample_id] = variant_com

if command == 'mosdepth':
    for sample_id in all_samples:
        tumor_id = os.path.join(input_dir, sample_id)
        output_prefix = os.path.join(output_dir, sample_id)

        mosdepth_com = [mosdepth, '--by', exonbed, output_prefix, tumor_id]
        submit_commands[sample_id] = conda_build + java_load + mosdepth_com

if __name__ == '__main__':
    # Submit jobs to cluster
    for sample_id, com in submit_commands.items():
        schedule_id = '{}_{}'.format(sample_id, command)
        arguments.schedule_job(command=com, name=schedule_id, python=python,
                               nodes=nodes, cores=cores, walltime=walltime)

