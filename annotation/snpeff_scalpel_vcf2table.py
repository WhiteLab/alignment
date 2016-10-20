#!/usr/bin/env python

import sys
import json
from utility.job_manager import job_manager


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['java'], config_data['tools']['snpsift'], config_data['params']['threads']


def convert_vcf(config_file, sample_pairs, suffix):
    (java, sift, th) = parse_config(config_file)
    cmd_list = []
    for pair in open(sample_pairs, 'r'):
        pair = pair.rstrip('\n').split('\t')
        pair = pair[0]
        in_vcf = pair + '/' + pair + suffix
        out_xls = pair + '/' + pair + '.indels.xls'
        cmd = java + ' -jar ' + sift + ' extractFields ' + in_vcf + \
              ' CHROM POS REF ALT "EFF[0].EFFECT" "EFF[0].FUNCLASS" "EFF[0].CODON" "EFF[0].AA" "EFF[0].AA_LEN" ' \
              '"EFF[0].GENE" "EFF[0].BIOTYPE" "EFF[0].CODING" MINCOV ALTCOV COVRATIO ID > ' + out_xls
        cmd_list.append(cmd)
    job_manager(cmd_list, th)
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='SNP annotation.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sp', '--sample_pairs', action='store', dest='sample_pairs', help='Sample tumor/normal pairs')
    parser.add_argument('-s', '--suffix', action='store', dest='suffix',
                        help='file suffix')


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample_pairs, suffix) = (
        inputs.config_file, inputs.sample_pairs,  inputs.suffix)
    convert_vcf(config_file, sample_pairs, suffix)
