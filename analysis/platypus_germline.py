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
        return config_data['tools']['platypus'], config_data['refs']['genome'], config_data['params']['threads']
    else:
        return config_data['tools']['platypus'], config_data['refs']['genome'], config_data['params']['threads'],\
               config_file['refs']['intervals'], config_file['params']['min_VAF_GL']


def platypus_germline(config_file, sample, log_dir, cflag):
    if cflag == 'y':
        (platypus, fasta, threads) = parse_config(config_file)
        platypus_cmd = platypus + " --nCPU=" + threads + " --refFile=" + fasta + " --bamFiles=" + sample \
                       + ".merged.final.bam -o " + sample + ".germline_calls/vcf > " + log_dir + sample \
                       + ".platypus.log 2>&1"
    else:
        (platypus, fasta, threads, region_file, minVAF) = parse_config(config_file)
        region_list = open(region_file, 'r')
        region_csv = region_list.read()
        regions = region_csv.replace('\n', ',')
        platypus_cmd = platypus + " --nCPU=" + threads + " --refFile=" + fasta + " --bamFiles=" + sample \
                       + ".merged.final.bam -o " + sample + ".germline_calls.vcf --minVarFreq=" + minVAF \
                       + "--regions=" + regions + " > " + log_dir + sample + ".platypus.log 2>&1"
    log(log_dir + sample + ".platypus.log", date_time() + platypus_cmd + "\n")
    f = 0
    try:
        f = subprocess.call(platypus_cmd, shell=True)
    except:
        log(log_dir + sample + ".platypus.log", 'platypus germline variant calling failed for sample ' + sample + '\n')
        exit(1)
    return f


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='novosort tool to sort BAM module.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/project name prefix')
    parser.add_argument('-l', '--log', action='store', dest='log_dir', help='LOG directory location')
    parser.add_argument('-f', '--flag', action='store', dest='cflag',
                        help='\'y\' if whole genome, \'n\' if custom capture to call only on-target regions')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample, log_dir, cflag) = (inputs.config_file, inputs.sample, inputs.log_dir, inputs.cflag)
    platypus_germline(config_file, sample, log_dir, cflag)
