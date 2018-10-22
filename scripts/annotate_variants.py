"""
Gregory Way 2017
scripts/annotate_variants.py

Use ANNOVAR to first convert a sample into annovar format and then annotate
"""

import os
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--merged', action='store_true',
                    help='use directory for merged VCFs')
parser.add_argument('-g', '--humanonly', action='store_true',
                    help='use humanonly directory for merged VCFs')
args = parser.parse_args()

merged = args.merged
humanonly = args.humanonly

if merged:
    vcf_file_dir = os.path.join('processed', 'gatk_merged_vcf')
    annovar_file_dir = os.path.join('results', 'annovar_merged_vcfs')
    annotated_file_dir = os.path.join('results', 'annotated_merged_vcfs')
else:
    vcf_file_dir = os.path.join('results', 'gatk_vcfs')
    annovar_file_dir = os.path.join('results', 'annovar_vcfs')
    annotated_file_dir = os.path.join('results', 'annotated_vcfs')

if humanonly:
    human_string = 'humanonly'
    vcf_file_dir = '{}_{}'.format(vcf_file_dir, human_string)
    annovar_file_dir = '{}_{}'.format(annovar_file_dir, human_string)
    annotated_file_dir = '{}_{}'.format(annotated_file_dir, human_string)

annovar_dir = os.path.join('modules', 'annovar')
humandb_dir = os.path.join(annovar_dir, 'humandb/')

convert_annovar = os.path.join(annovar_dir, 'convert2annovar.pl')
table_annovar = os.path.join(annovar_dir, 'table_annovar.pl')

conv_com = 'perl {} -format vcf4 -filter pass'.format(convert_annovar)
anno_com = 'perl {} {} -buildver hg19'.format(table_annovar, humandb_dir)

# Convert to annovar format
for vcf_file in os.listdir(vcf_file_dir):
    if '.idx' not in vcf_file:
        base_name = vcf_file.split('.')[0]
        full_vcf_file = os.path.join(vcf_file_dir, vcf_file)
        output_vcf_file = os.path.join(annovar_file_dir,
                                       '{}.annovar.vcf'.format(base_name))
        if not os.path.isfile(output_vcf_file):
            file_command = '{} {} > {}'.format(conv_com, full_vcf_file,
                                               output_vcf_file)
            subprocess.call(file_command, shell=True)

# Annotate annovar formatted files with given databases
for annovar_file in os.listdir(annovar_file_dir):
    base_name = annovar_file.split('.')[0]
    full_annov_file = os.path.join(annovar_file_dir, annovar_file)
    annotated_vcf_file = os.path.join(annotated_file_dir,
                                      '{}.annotated'.format(base_name))
    if not os.path.isfile(annotated_vcf_file):
        file_command = 'perl {} {} modules/annovar/humandb -buildver hg19 ' \
                       '-out {} -verbose -otherinfo -remove -protocol ' \
                       'refGene,cosmic70,gnomad_exome,dbnsfp30a ' \
                       '-operation g,f,f,f -nastring . -csvout ' \
                       '-polish'.format(table_annovar, full_annov_file,
                                        annotated_vcf_file)
        subprocess.call(file_command, shell=True)
