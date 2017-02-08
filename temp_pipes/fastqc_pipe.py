#!/usr/bin/env python

import sys
import json
import os
import glob
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.job_manager import job_manager


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['refs']['cont'], config_data['refs']['obj'], config_data['tools']['fastqc'],\
           config_data['params']['threads']


def fastqc_pipe(flist, config_file):
    (cont, obj, fastqc_tool, threads) = parse_config(config_file)
    src_cmd = '. /home/ubuntu/.novarc;'
    job_list = []
    for fq in open(flist):
        fq = fq.rstrip('\n')
        root = os.path.basename(fq).replace('_sequence.txt.gz', '')
        parts = fq.split('/')
        dl_cmd = src_cmd + 'swift download ' + cont + ' ' + fq
        outdir = obj + '/' + parts[1] + '/QC/'
        logdir = obj + '/' + parts[1] + '/LOGS/'
        setup_cmd = 'mkdir -p ' + outdir + ' ' + logdir + ';'
        logfile = logdir + root + '.fastqc.log'
        fastqc_cmd = fastqc_tool + ' -o ' + outdir + ' ' + fq + '2> ' + logfile + ';'
        up_cmd = src_cmd + 'swift upload ' + cont + ' ' + logfile + ';'
        up_list = glob.glob(outdir + root + '*')
        for fn in up_list:
            up_cmd += src_cmd + 'swift upload ' + cont + ' ' + fn + ';'
        cleanup = 'rm ' + fq + ';'
        final_cmd = dl_cmd + setup_cmd + fastqc_cmd + up_cmd + cleanup
        job_list.append(final_cmd)
    job_manager(job_list, threads)


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Temp pipeline to run FastQC over sequnecing files and upload results')
    parser.add_argument('-f', '--fastq_list', action='store', dest='flist', help='List of fastq files')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()

    flist = inputs.flist
    config_file = inputs.config_file
    fastqc_pipe(flist, config_file)


if __name__ == "__main__":
    main()
