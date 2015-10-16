#!/usr/bin/env python
__author__ = 'Miguel'
import sys
import json

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from job_manager import job_manager


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['samtools'], config_data['tools']['varscan'], config_data['refs']['fai'], \
           config_data['refs']['fa_ordered'], config_data['params']['threads']


def varscan_germline(config_file, sample):
    (samtools, varscan, region, fasta, th) = parse_config(config_file)
    rf = open(region, 'r')
    cmd_list = []
    for line in rf:
        chrom = line.split('\t')
        cmd = samtools + ' mpileup -r ' + chrom[
            0] + ' -B -f ' + fasta + ' ' + sample + '.merged.final.bam | java -jar ' + varscan + ' mpileup2cns'
        cmd_list.append(cmd)
    rf.close()
    proc = int(th) - 2
    job_manager(cmd_list, str(proc))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Germline variant caller.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample prefix')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample) = (inputs.config_file, inputs.sample)
    varscan_germline(config_file, sample)
