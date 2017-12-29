#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
import subprocess
from utility.log import log
import json
import os


def parse_config(config_file, cflag):
    config_data = json.loads(open(config_file, 'r').read())
    if cflag == 'y':
        return config_data['tools']['platypus'], config_data['refs']['fa_ordered'], config_data['params']['threads'], \
               config_data['refs']['project_dir'], config_data['refs']['project'], config_data['refs']['align_dir']
    else:
        return config_data['tools']['platypus'], config_data['refs']['fa_ordered'], config_data['params']['threads'],\
               config_data['refs']['intervals_0base'], config_data['params']['min_VAF_GL'], \
               config_data['tools']['samtools'], config_data['refs']['project_dir'], config_data['refs']['project'], \
               config_data['refs']['align_dir']


def platypus_germline(config_file, sample, log_dir, cflag):

    loc = log_dir + sample + ".platypus.log"
    # here for safety as python is confusing about whether variables exist outside of if-else statements or not
    platypus_cmd = ''
    if cflag == 'y':
        (platypus, fasta, threads, project_dir, project, align) = parse_config(config_file, cflag)
        bam = project_dir + project + '/' + align + '/' + sample + '/BAM/' + sample + '.merged.final.bam'
        platypus_cmd = "python2.7 " + platypus + " callVariants --nCPU=" + threads + " --refFile=" + fasta \
                       + " --bamFiles=" + bam + " -o " + sample + ".germline_calls.vcf --logFileName=" \
                       + log_dir + sample + ".platypus.log" + " >> " + loc + " 2>&1"
    else:
        (platypus, fasta, threads, region_file, minVAF, samtools, project_dir, project, align) \
            = parse_config(config_file, cflag)

        bam = project_dir + project + '/' + align + '/' + sample + '/BAM/' + sample + '.merged.final.bam'
        if not (os.path.isfile(bam + '.bai') or os.path.isfile(bam[:-1] + 'i')):
            log(loc, date_time() + bam + ' not indexed.  Indexing\n')
            cmd = samtools + ' index ' + bam
            log(loc, date_time() + cmd + '\n')
            subprocess.call(cmd, shell=True)
        platypus_cmd = "python2.7 " + platypus + " callVariants --nCPU=" + threads + " --refFile=" + fasta \
                       + " --bamFiles=" + bam + " --filterDuplicates=0 -o " + sample \
                       + ".germline_calls.vcf --minVarFreq=" + minVAF + " --regions=" + region_file \
                       + " --logFileName=" + loc + " >> " + loc + " 2>&1"
    log(loc, date_time() + platypus_cmd + "\n")
    f = 0
    try:
        f = subprocess.call(platypus_cmd, shell=True)
    except:
        log(loc, 'platypus germline variant calling failed for sample ' + sample + '\n')
        return f

    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Germline calling using Platypus.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-s', '--sample', action='store', dest='sample', help='Normal sample ID')
    parser.add_argument('-l', '--log', action='store', dest='log_dir', help='LOG directory location')
    parser.add_argument('-f', '--flag', action='store', dest='cflag',
                        help='\'y\' if whole genome, \'n\' if custom capture to call only on-target regions')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample, log_dir, cflag) = (inputs.config_file, inputs.sample, inputs.log_dir, inputs.cflag)
    platypus_germline(config_file, sample, log_dir, cflag)
