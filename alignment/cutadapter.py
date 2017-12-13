#!/usr/bin/env python3
import os
import json
import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from subprocess import call
from subprocess import Popen
from utility.log import log
from alignment.fastqc import fastqc
import time

def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    (cutadapt_tool, fastqc_tool, qual, mqual, r1_adapt, r2_adapt, minlen, r1trim, r2trim, aflag, ntrim, threads) = (
        config_data['tools']['cutadapt'], config_data['tools']['fastqc'], config_data['params']['qual'],
        config_data['params']['mqual'], config_data['params']['r1adapt'], config_data['params']['r2adapt'],
        config_data['params']['minlen'], config_data['params']['r1trim'], config_data['params']['r2trim'],
        config_data['params']['aflag'], config_data['params']['ntrim'], config_data['params']['threads'])
    return cutadapt_tool, fastqc_tool, qual, mqual, r1_adapt, r2_adapt, minlen, r1trim, r2trim, aflag, ntrim, threads


def cutadapter(sample, end1, end2, config_file):
    # casual logging - look for a LOGS directory, otherwise assume current dir
    log_dir = './'
    if os.path.isdir('LOGS'):
        log_dir = 'LOGS/'
    loc = log_dir + sample + '.cutadapt.log'
    temp1 = end1 + '.temp.gz'
    temp2 = end2 + '.temp.gz'
    (cutadapt_tool, fastqc_tool, qual, mqual, r1_adapt, r2_adapt, minlen, r1trim, r2trim, aflag, ntrim, threads) = \
        parse_config(config_file)
    aflag2 = aflag.upper()
    # cut_th = str(int(threads) - 2)
    cut_th = "3"
    cutadapt_cmd = cutadapt_tool + ' -j ' + cut_th + ' -m ' + minlen + ' --quality-base=' + qual + ' -q ' + mqual \
                   + ' -' + aflag + ' ' + r1_adapt + ' -' + aflag2 + ' ' + r2_adapt + ' -u ' + r1trim + ' -U ' \
                   + r2trim + ' -n ' + ntrim + ' -o ' + temp1 + ' -p ' + temp2 + ' ' + end1 + ' ' + end2 + ' >> ' \
                   + loc + ' 2>> ' + loc
    # cutadapt 1.15 can now quality trim without having to fake an adapter
    if r1_adapt[0] == 'Z':
        cutadapt_cmd = cutadapt_tool + ' -j ' + cut_th + ' -m ' + minlen + ' --quality-base=' + qual + ' -q ' + mqual \
                       + ' -u ' + r1trim + ' -U ' + r2trim + ' -n ' + ntrim + ' -o ' + temp1 + ' -p ' + temp2 + ' ' \
                       + end1 + ' ' + end2 + ' >> ' + loc + ' 2>> ' + loc
    log(loc, date_time() + cutadapt_cmd + '\n')
    f = Popen(cutadapt_cmd, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
    time.sleep(20)

    if str(f.poll()) == '1':
        log(loc, date_time() + 'cutadapt returned an error.  Check your inputs and try again!\n')
        exit(1)
    else:
        log(loc, date_time() + 'cutadapt running ok.  Starting FastQC\n')
    # will run fastqc while cutadapt is running - assuming a vm of at least 4 cores
    check = fastqc(fastqc_tool, sample, end1, end2, '2')
    if check != 0:
        log(loc, date_time() + 'FastQC failed for sample! ' + sample + '\n')
        exit(1)
    else:
        log(loc, date_time() + 'FastQC completed running, checking on cutadapt\n')
    # poll cutadapt until complete
    while f.poll() is None:
        time.sleep(10)
    log(loc, date_time() + 'Cutadapt completed running, renaming files\n')
    rn_fq = 'mv ' + temp1 + ' ' + end1 + ';mv ' + temp2 + ' ' + end2
    call(rn_fq, shell=True)
    return 0


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description='cutadapt module.  Removes Set minimum base quality')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/location name prefix')
    parser.add_argument('-f1', '--file1', action='store', dest='end1', help='First of paired-end fastq file')
    parser.add_argument('-f2', '--file2', action='store', dest='end2', help='Second of paired-end fastq file')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample, end1, end2, config_file) = (inputs.sample, inputs.end1, inputs.end2, inputs.config_file)
    cutadapter(sample, end1, end2, config_file)
