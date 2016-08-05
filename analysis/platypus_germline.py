#!/usr/bin/env python
import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
import subprocess
from log import log
import json


def parse_config(config_file, cflag):
    config_data = json.loads(open(config_file, 'r').read())
    if cflag == 'y':
        return config_data['tools']['platypus'], config_data['refs']['fa_ordered'], config_data['params']['threads']
    else:
        return config_data['tools']['platypus'], config_data['refs']['fa_ordered'], config_data['params']['threads'],\
               config_data['refs']['intervals_0base'], config_data['params']['min_VAF_GL']


def platypus_germline(config_file, samp_list, log_dir, cflag, ref_mnt):
    samp_fh = open(samp_list, 'r')
    # only want to run germline once if norm being used multiple times as a comparison
    samp_flag = {}
    for line in samp_fh:
        info = line.rstrip('\n').split('\t')
        sample = info[0]
        if len(info) > 1:
            sample = info[2]
        loc = log_dir + sample + ".platypus.log"
        if sample not in samp_flag:

            if cflag == 'y':
                (platypus, fasta, threads) = parse_config(config_file, cflag)
                fasta = ref_mnt + '/' + fasta
                platypus_cmd = platypus + " callVariants --nCPU=" + threads + " --refFile=" + fasta + " --bamFiles=" + sample \
                               + ".merged.final.bam -o " + sample + ".germline_calls.vcf --logFileName=" + log_dir + sample \
                               + ".platypus.log" + " >> " + log_dir + sample + ".platypus.log 2>&1"
            else:
                (platypus, fasta, threads, region_file, minVAF) = parse_config(config_file, cflag)
                fasta = ref_mnt + '/' + fasta
                regions = ref_mnt + '/' + region_file
                platypus_cmd = platypus + " callVariants --nCPU=" + threads + " --refFile=" + fasta + " --bamFiles=" + sample \
                               + ".merged.final.bam --filterDuplicates=0 -o " + sample + ".germline_calls.vcf --minVarFreq="\
                               + minVAF + " --regions=" + regions + " --logFileName=" + log_dir + sample \
                               + ".platypus.log >> " + log_dir + sample + ".platypus.log 2>&1"
            log(log_dir + sample + ".platypus.log", date_time() + platypus_cmd + "\n")
            f = 0
            try:
                f = subprocess.call(platypus_cmd, shell=True)
                samp_flag[samp_fh] = 1
            except:
                log(loc, 'platypus germline variant calling failed for sample ' + sample
                    + '\n')
                return f
        else:
            log(loc, 'platypus already run on ' + sample + ' during this session, skipping!\n')
        return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Germline calling using Platypus.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sp', '--sample_pairs', action='store', dest='sample_pairs', help='Sample tumor/normal pairs,'
                                                                                           ' or sample list')
    parser.add_argument('-l', '--log', action='store', dest='log_dir', help='LOG directory location')
    parser.add_argument('-f', '--flag', action='store', dest='cflag',
                        help='\'y\' if whole genome, \'n\' if custom capture to call only on-target regions')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference drive path - i.e. /mnt/cinder/REFS_XXXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample, log_dir, cflag, ref_mnt) = (inputs.config_file, inputs.sample_pairs, inputs.log_dir, inputs.cflag,
                                                      inputs.ref_mnt)
    platypus_germline(config_file, sample, log_dir, cflag, ref_mnt)
