# Gregory Way
# Modified from @jietan of bitbucket.org/greenelab

import sys
import argparse
from subprocess import call


class PBSJob:
    """
    Job for working with Torque at Dartmouth
    If you use job-arrays, the dartmouth-defined limit is 100, for more
    you have to submit jobs individually
    """
    def __init__(self, name='', queue='default', nodes=1, ppn=1,
                 walltime='01:00:00', mail='a', addr='gregway@mail.med.upenn.edu',
                 cwd=True, command=None, array=None):
        self.name = name
        self.queue = queue
        self.nodes = nodes
        self.ppn = ppn
        self.walltime = walltime
        self.mail = mail
        self.array = array
        if addr is None:
            sys.stderr.write('An email address is REQUIRED')
            return None
        self.addr = addr
        self.cwd = cwd
        if command is None:
            sys.stderr.write('A command is REQUIRED')
            return None
        self.command = command

    def set_command(self, command):
        self.command = command

    def set_name_command(self, name, command):
        self.name = name
        self.command = command

    def write(self, filename):
        with open(filename, 'w') as fh:
            fh.write('#!/bin/bash -l\n')
            fh.write('#PBS -N ' + str(self.name) + '\n')
            fh.write('#PBS -q ' + str(self.queue) + '\n')
            fh.write('#PBS -l nodes=' + str(self.nodes) + ':ppn=' + str(self.ppn) + '\n')
            fh.write('#PBS -l walltime=' + str(self.walltime) + '\n')
            fh.write('#PBS -m ' + str(self.mail) + '\n')
            fh.write('#PBS -M ' + str(self.addr) + '\n')
            if self.array is not None:
                fh.write('#PBS -t 1-' + str(self.array) + '\n')
            fh.write('export TERM=xterm\n')
            if self.cwd:
                fh.write('cd $PBS_O_WORKDIR\n')
            fh.write(str(self.command) + '\n')
            fh.write('exit 0\n')

# Command arguments to submit a job from command line
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--command', help='what command to submit to cluster')
parser.add_argument('-n', '--name', help='the name of the command to call')
parser.add_argument('-w', '--walltime', help='the amount of time alloted to the script')
parser.add_argument('-f', '--filename', help='what file name to save the PBS script')
parser.add_argument('-o', '--nodes', default=1, help='the number of nodes to allocate')
parser.add_argument('-r', '--cores', default=1, help='the number of cores to allocate per node')

args = parser.parse_args()
COMMAND = args.command
NAME = args.name
WALLTIME = args.walltime
FILENAME = args.filename
NODES = args.nodes
CORES = args.cores

# Generate the script and run the command
job = PBSJob(command=COMMAND, name=NAME, walltime=WALLTIME, nodes=NODES, ppn=CORES)
job.write(FILENAME)
call(['qsub', FILENAME])



