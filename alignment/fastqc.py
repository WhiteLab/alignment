#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
import os
from utility.date_time import date_time
from subprocess import Popen
from subprocess import call
from utility.log import log


def fastqc(fastqc_tool, sample, end1, end2, t):
    # casual logging - look for a LOGS directory, otherwise assume current dir
    log_dir = './'
    if os.path.isdir('LOGS'):
        log_dir = 'LOGS/'
    loc = log_dir + sample + '.fastqc.log'
    fastqc_cmd = fastqc_tool + ' --extract -t ' + t + ' -o QC/ ' + end1 + ' ' + end2
    log(loc, date_time() + fastqc_cmd + "\n")
    f = Popen(fastqc_cmd, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    # check after a minute whether the process is still good - shouldn't take too long to ascertain whether phred score
    #  didn't fit
    call('sleep 20s', shell=True)

    if str(f.poll()) == '1':
        log(loc, date_time() + 'fastqc returned an error.  Check your inputs and try again!\n')
        exit(1)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='fastqc module.  Provides quality stats for fastq file and is independent of alignment.')
    parser.add_argument('-f', '--fastqc', action='store', dest='fastqc_tool', help='Location of fastqc tool.')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/location name prefix')
    parser.add_argument('-f1', '--file1', action='store', dest='end1', help='First of paired-end fastq file')
    parser.add_argument('-f2', '--file2', action='store', dest='end2', help='Second of paired-end fastq file')
    parser.add_argument('-t', '--threads', action='store', dest='t', help='Number of threads')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (fastqc_tool, sample, end1, end2, t) = (inputs.fastqc_tool, inputs.sample, inputs.end1, inputs.end2, inputs.t)
    fastqc(fastqc_tool, sample, end1, end2, t)
