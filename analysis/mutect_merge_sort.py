#!/usr/bin/python
import sys
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


def mutect_merge_sort(config_file, sample_pairs, ref_mnt):
    # use fasta index to get sort order for file output
    (fai) = parse_config(config_file)
    fai_list = []
    fai_fh = open(ref_mnt + '/' + fai, 'r')
    for line in fai_fh:
        line = line.rstrip('\n')
        data = line.split('\t')
        fai_list.append(data[0])
    fai_fh.close()
    # output files should be in directory named after sample-pairs
    sp_fh = open(sample_pairs, 'r')
    for line in sp_fh:
        line = line.rstrip('\n')
        sp = line.split('\t')
        dir_list = os.listdir(sp[0])
        suffix_dict = {}
        for fn in dir_list:
            parts = fn.split('.')
            if parts[2] == 'out' or parts[2] == 'vcf':
                suffix = '.'.join(parts[2:])
                if suffix not in suffix_dict:
                    suffix_dict[suffix] = []
                suffix_dict[suffix].append(fn)
        merge_sort(suffix_dict, sp[0], fai_list)
    sys.stderr.write(date_time() + 'File merging completed\n')
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Merge and sort output from mutect variant caller.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sp', '--sample_pairs', action='store', dest='sample_pairs', help='Sample tumor/normal pairs')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference drive path - i.e. /mnt/cinder/REFS_XXXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample_pairs, ref_mnt) = (inputs.config_file, inputs.sample_pairs, inputs.ref_mnt)
    mutect_merge_sort(config_file, sample_pairs, ref_mnt)
