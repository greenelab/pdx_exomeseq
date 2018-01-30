"""
Gregory Way 2017
scripts/util/arguments.py

Function to facilitate argument creation - Import only
"""

import os
import re
import argparse
import subprocess


def get_args():
    parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(title='subparsers', dest='command')

    # Create universal command arguments
    parser.add_argument('-w', '--walltime', default='04:00:00',
                        help='the amount of time alloted to the script')
    parser.add_argument('-n', '--nodes', default=1,
                        help='the number of nodes to allocate')
    parser.add_argument('-r', '--cores', default=4,
                        help='the number of cores to allocate per node')
    parser.add_argument('-i', '--input_directory', default=None,
                        help='the directory of the files to process')
    parser.add_argument('-o', '--output_directory', default='.',
                        help='the directory to save the output files')
    parser.add_argument('-y', '--config_yaml',
                        default='discovery_variables.yml',
                        help='Configuration variables for input')
    parser.add_argument('-g', '--genome', default='hg',
                        help='name of the reference genome')

    # Create subcommands for specific pipeline functionality
    parser_fastqc = subparsers.add_parser('fastqc', parents=[parser])
    parser_fastqc.set_defaults(func=get_fastqc, which='fastqc')

    parser_multiqc = subparsers.add_parser('multiqc', parents=[parser])
    parser_multiqc.add_argument('--html_file')
    parser_multiqc.set_defaults(func=get_fastqc, which='multiqc')

    parser_trimgalore = subparsers.add_parser('trimgalore', parents=[parser])
    parser_trimgalore.add_argument('--fastqc_results_dir')
    parser_trimgalore.set_defaults(func=get_trimgalore, which='trimgalore')

    parser_bwa = subparsers.add_parser('bwa', parents=[parser])
    parser_bwa.set_defaults(func=get_bwa, which='bwa')

    parser_samtools = subparsers.add_parser('samtools', parents=[parser])
    parser_samtools.add_argument('--sub_command')
    parser_samtools.set_defaults(func=get_samtools, which='samtools')

    parser_disambiguate = subparsers.add_parser('disambiguate',
                                                parents=[parser])
    parser_disambiguate.add_argument('--human_dir')
    parser_disambiguate.add_argument('--mouse_dir')
    parser_disambiguate.set_defaults(func=get_disambiguate,
                                     which='disambiguate')

    parser_variant = subparsers.add_parser('variant', parents=[parser])
    parser_variant.add_argument('--sub_command')
    parser_variant.add_argument('--threads', default='6')
    parser_variant.add_argument('--min_confidence', default='20')
    parser_variant.add_argument('--min_base_quality_score', default='20')
    parser_variant.set_defaults(func=get_variant, which='variant')

    parser_mosdepth = subparsers.add_parser('mosdepth', parents=[parser])
    parser_mosdepth.set_defaults(func=get_mosdepth, which='mosdepth')

    args = parser.parse_args()
    return args


def schedule_job(command, name, python, nodes=1, cores=4, walltime='04:00:00'):
    output_com = [python, os.path.join('util', 'schedule.py'),
                  '--command', ' '.join(str(x) for x in command),
                  '--name', name,
                  '--walltime', walltime,
                  '--nodes', str(nodes),
                  '--cores', str(cores),
                  '--filename', os.path.join('logs', '{}.pbs'.format(name))]
    return subprocess.call(output_com)


def get_fastqc(args):
    wes_files = []
    for path, subdirs, files in os.walk(args.input_directory):
        for name in files:
            if re.match(r'(fasta|fq|fastq)\.gz$', name):
                full_name = os.path.join(path, name)
                wes_files.append(full_name)
    return wes_files


def get_trimgalore(args):
    wes_files = []
    for path, subdirs, files in os.walk(args.input_directory):
        for name in files:
            if 'fastq.gz' in name:
                full_name = os.path.join(path, name)
                wes_files.append(full_name)
    paired_reads = []
    for name in wes_files:
        if '_R1_' in name:
            read_1 = name
            read_2 = name.replace('_R1_', '_R2_')
            paired_reads.append([read_1, read_2])
    return paired_reads


def get_bwa(args):
    all_files = []
    for path, subdirs, files in os.walk(args.input_directory):
        for name in files:
            if 'fq.gz' in name:
                all_files.append(name)

    paired_reads = []
    for name in all_files:
        if '_R1_' in name:
            read_1 = name
            read_2 = name.replace('_R1_', '_R2_')
            read_2 = read_2.replace('val_1', 'val_2')
            paired_reads.append([read_1, read_2])
    return paired_reads


def get_samtools(args):
    check_suffix = ['fq.gz.sam', '.sam_sorted.bam', '.bam_sorted_fixmate.bam',
                    'disambiguatedSpeciesA.bam', 'merged.bam']
    sam_files = []
    for path, subdirs, files in os.walk(args.input_directory):
        for name in files:
            if (any([x in name for x in check_suffix])) and ('.bai' not in name):
                sam_files.append(name)
    return sam_files


def get_disambiguate(args):
    sample_files = []
    for path, subdirs, files in os.walk(args.human_dir):
        for name in files:
            sample_id = name.split('.')[0]  # This extracts the sample ID
            sample_files.append(sample_id)
    return sample_files


def get_variant(args):
    check_suffix = ['bam_rmdup.bam', 'merged.bam']
    bam_files = []
    for path, subdirs, files in os.walk(args.input_directory):
        for name in files:
            if (any([x in name for x in check_suffix])) and ('.bai' not in name):
                bam_files.append(name)
    return bam_files


def get_mosdepth(args):
    bam_files = []
    for path, subdirs, files in os.walk(args.input_directory):
        for name in files:
            if '.bai' not in name:
                bam_files.append(name)
    return bam_files
