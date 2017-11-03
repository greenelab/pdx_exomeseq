import os
import subprocess


def schedule_job(command, name, step, nodes=1, cores=4, walltime='04:00:00'):
    output_com = ['python', os.path.join('util', 'schedule.py'),
                  '--command', command,
                  '--name', name,
                  '--walltime', walltime,
                  '--nodes', str(nodes),
                  '--cores', str(cores),
                  '--filename', 'logs/' + step + '_' + name + '.pbs']
    return output_com


all_files = []
for path, subdirs, files in os.walk('data/trimmed'):
    for name in files:
        if 'trimmed.fq.gz' in name:
            all_files.append(name)

paired_reads = []
for name in all_files:
    if '_R1_' in name:
        read_1 = name
        read_2 = name.replace('_R1_', '_R2_') 
        paired_reads.append(['data/trimmed/' + read_1, 
                             'data/trimmed/' + read_2])

#for pair_combo in paired_reads:
for sample_1, sample_2 in paired_reads:
    
    name = os.path.basename(sample_1)

    # Load constants
    base_dir = '/lorax/sanchezlab/shared/pdx_exomeseq/'
    bwa = base_dir + 'modules/bwa-0.7.5a/bwa'
    hg_ref_file = base_dir + 'reference/ucsc.hg19.fasta'
    mm_ref_file = base_dir + 'reference/mm9.fa'
    sample_basename = base_dir + 'data/bwa/' + name.replace('_R1_', '_')

    # Output files
    sample_1_hg_sai = sample_basename + '_hg_1.sai'
    sample_2_hg_sai = sample_basename + '_hg_2.sai'
    sample_1_mm_sai = sample_basename + '_mm_1.sai'
    sample_2_mm_sai = sample_basename + '_mm_2.sai'

    sample_hg_sam = sample_basename + '_hg_aln.sam'
    sample_mm_sam = sample_basename + '_mm_aln.sam'

    sample_hg_bam = sample_basename + '_hg.bam'
    sample_mm_bam = sample_basename + '_mm.bam'

    disambig_dir = base_dir + 'data/disambiguated/'

    # Generate the command calls
    bwa_1_hg_com = bwa + ' aln ' + hg_ref_file + ' ' + sample_1 + ' > ' + sample_1_hg_sai
    bwa_2_hg_com = bwa + ' aln ' + hg_ref_file + ' ' + sample_2 + ' > ' + sample_2_hg_sai
    bwa_1_mm_com = bwa + ' aln ' + mm_ref_file + ' ' + sample_1 + ' > ' + sample_1_mm_sai
    bwa_2_mm_com = bwa + ' aln ' + mm_ref_file + ' ' + sample_2 + ' > ' + sample_2_mm_sai

    bwa_paired_hg_com = bwa + ' sampe ' + hg_ref_file + ' ' + \
                    sample_1_hg_sai + ' ' + sample_2_hg_sai + ' ' + \
                    sample_1 + ' ' + sample_2 + ' > ' + sample_hg_sam
    bwa_paired_mm_com = bwa + ' sampe ' + mm_ref_file + ' ' + \
                    sample_1_mm_sai + ' ' + sample_2_mm_sai + ' ' + \
                    sample_1 + ' ' + sample_2 + ' > ' + sample_mm_sam

    samtools_sorted_bam_hg_com = 'samtools view -bS ' + sample_hg_sam + ' > ' + sample_hg_bam
    samtools_sorted_bam_mm_com = 'samtools view -bS ' + sample_mm_sam + ' > ' + sample_mm_bam

    disambiguate_com = '~/.conda/envs/pdx-exomeseq/bin/ngs_disambiguate ' + \
                   sample_hg_bam + ' ' + sample_mm_bam + \
                   ' --prefix test --aligner bwa --output-dir ' + disambig_dir

    # Make the commands in sequence
    bwa_1_hg_sched = schedule_job(bwa_1_hg_com, name, 'bwa_hg1_step')
    subprocess.call(bwa_1_hg_sched)

    bwa_2_hg_sched = schedule_job(bwa_2_hg_com, name, 'bwa_hg2_step')
    subprocess.call(bwa_2_hg_sched)

    bwa_1_mm_sched = schedule_job(bwa_1_mm_com, name, 'bwa_mm1_step')
    subprocess.call(bwa_1_mm_sched)

    bwa_2_mm_sched = schedule_job(bwa_2_mm_com, name, 'bwa_mm2_step')
    subprocess.call(bwa_2_mm_sched)

    #bwa_paired_hg_sched = schedule_job(bwa_paired_hg_com, name, 'bwa_hg_sampe')
    #subprocess.call(bwa_paired_hg_sched)

    #bwa_paired_mm_sched = schedule_job(bwa_paired_mm_com, name, 'bwa_mm_sampe')
    #subprocess.call(bwa_paired_mm_sched)

    #samtools_sorted_hg_sched = schedule_job(samtools_sorted_bam_hg_com, name, 'samtools_hg_sort')
    #subprocess.call(samtools_sorted_hg_sched)

    #samtools_sorted_mm_sched = schedule_job(samtools_sorted_bam_mm_com, name, 'samtools_mm_sort')
    #subprocess.call(samtools_sorted_mm_sched)

    #disambiguate_sched = schedule_job(disambiguate_com, name, 'disambiguate')
    #subprocess.call(disambiguate_sched)


