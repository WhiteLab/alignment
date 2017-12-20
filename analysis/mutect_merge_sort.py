#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
import json
import os


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['refs']['fai']


def merge_sort(suffix_dict, pair, fai_list):
    sys.stderr.write(date_time() + 'Started creating merged out file for ' + pair + '\n')
    if 'out' in suffix_dict:
        out = {}
        head = ''
        for fn in suffix_dict['out']:
            fh = open(pair + '/' + fn, 'r')
            head = next(fh)
            head += next(fh)
            for line in fh:
                data = line.split('\t')
                data[1] = int(data[1])
                if data[0] not in out:
                    out[data[0]] = {}
                if data[1] not in out[data[0]]:
                    out[data[0]][data[1]] = []
                out[data[0]][data[1]].append(line)
            fh.close()
        fh = open(pair + '.out', 'w')
        fh.write(head)
        for chrom in fai_list:
            if chrom in out:
                for pos in sorted(out[chrom]):
                    for line in out[chrom][pos]:
                        fh.write(line)
        fh.close()
        sys.stderr.write(date_time() + 'Started creating merged out.keep file for ' + pair + '\n')
    if 'out.keep' in suffix_dict:
        out = {}
        head = ''
        for fn in suffix_dict['out.keep']:
            fh = open(pair + '/' + fn, 'r')
            head = ''
            head = next(fh)
            head += next(fh)
            for line in fh:
                data = line.split('\t')
                data[1] = int(data[1])
                if data[0] not in out:
                    out[data[0]] = {}
                if data[1] not in out[data[0]]:
                    out[data[0]][data[1]] = []
                out[data[0]][data[1]].append(line)
            fh.close()
        fh = open(pair + '.out.keep', 'w')
        fh.write(head)
        for chrom in fai_list:
            if chrom in out:
                for pos in sorted(out[chrom]):
                    for line in out[chrom][pos]:
                        fh.write(line)
        fh.close()
    sys.stderr.write(date_time() + 'Started creating merged vcf file for ' + pair + '\n')
    if 'vcf' in suffix_dict:
        out = {}
        head = ''
        for fn in suffix_dict['vcf']:
            fh = open(pair + '/' + fn, 'r')
            head = ''
            for line in fh:
                head += line
                if line[0:2] != '##':
                    break
            for line in fh:
                data = line.split('\t')
                data[1] = int(data[1])
                if data[0] not in out:
                    out[data[0]] = {}
                if data[1] not in out[data[0]]:
                    out[data[0]][data[1]] = []
                out[data[0]][data[1]].append(line)
            fh.close()
        fh = open(pair + '.vcf', 'w')
        fh.write(head)
        for chrom in fai_list:
            if chrom in out:
                for pos in sorted(out[chrom]):
                    for line in out[chrom][pos]:
                        fh.write(line)
        fh.close()
    sys.stderr.write(date_time() + 'Started creating merged vcf.keep file for ' + pair + '\n')
    if 'vcf.keep' in suffix_dict:
        out = {}
        head = ''
        for fn in suffix_dict['vcf.keep']:
            fh = open(pair + '/' + fn, 'r')
            head = ''
            for line in fh:
                head += line
                if line[0:2] != '##':
                    break
            for line in fh:
                data = line.split('\t')
                data[1] = int(data[1])
                if data[0] not in out:
                    out[data[0]] = {}
                if data[1] not in out[data[0]]:
                    out[data[0]][data[1]] = []
                out[data[0]][data[1]].append(line)
            fh.close()
        fh = open(pair + '.vcf.keep', 'w')
        fh.write(head)
        for chrom in fai_list:
            if chrom in out:
                for pos in sorted(out[chrom]):
                    for line in out[chrom][pos]:
                        fh.write(line)
        fh.close()

    sys.stderr.write(date_time() + 'File merging complete for ' + pair + '\n')


def mutect_merge_sort(config_file, sample_pair):
    # use fasta index to get sort order for file output
    (fai) = parse_config(config_file)
    fai_list = []
    fai_fh = open(fai, 'r')
    for line in fai_fh:
        line = line.rstrip('\n')
        data = line.split('\t')
        fai_list.append(data[0])
    fai_fh.close()
    # output files should be in directory named after sample-pairs

    dir_list = os.listdir('./')
    suffix_dict = {}
    for fn in dir_list:
        parts = fn.split('.')
        if len(parts) >= 3:
            if parts[2] == 'out' or parts[2] == 'vcf':
                suffix = '.'.join(parts[2:])
                if suffix not in suffix_dict:
                    suffix_dict[suffix] = []
                suffix_dict[suffix].append(fn)
    merge_sort(suffix_dict, sample_pair, fai_list)
    sys.stderr.write(date_time() + 'File merging completed\n')
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Merge and sort output from mutect variant caller.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-s', '--sample_pair', action='store', dest='sample_pair', help='Sample tumor/normal pairs')


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample_pair) = (inputs.config_file, inputs.sample_pair)
    mutect_merge_sort(config_file, sample_pair)
