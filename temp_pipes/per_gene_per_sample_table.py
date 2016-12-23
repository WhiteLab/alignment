#!/usr/bin/env python

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.date_time import date_time


def process_bed(bed):
    gene_list = []
    gene_dict = {}
    master_dict = {}
    for entry in open(bed):
        info = entry.rstrip('\n').split('\t')
        parts = info[-1].split('_')
        if parts[0] not in gene_dict:
            gene_list.append(parts[0])
            gene_dict[parts[0]] = {}
            master_dict[parts[0]] = {}
        gene_dict[parts[0]] = {'len': (int(info[1]) - int(info[2])), 'tot_cov': 0}
    return gene_list, gene_dict, master_dict


def calc_coverage(sample_list, suffix, bed):
    slist = []
    sys.stderr.write(date_time() + 'Processing bed file ' + bed + '\n')
    (gene_list, gene_dict, master_dict) = process_bed(bed)
    for sample in open(sample_list):
        sys.stderr.write(date_time() + 'Processing sample ' + sample)
        slist.append(sample)
        sample = sample.rstrip('\n')
        cur = sample + suffix
        temp_dict = gene_dict
        for entry in open(cur):
            info = entry.rstrip('\n').split('\t')
            if info[0] == 'all':
                break
            parts = info[3].split('_')
            temp_dict[parts[0]]['tot_cov'] += (int(info[4]) * int(info[5]))
        for gene in gene_list:
            master_dict[gene][sample] = (float(temp_dict[gene]['tot_cov'])/temp_dict[gene]['len'])
    sys.stderr.write(date_time() + 'Outputting results\n')
    sys.stdout.write('Gene/Sample')
    print '\t'.join(slist)
    for gene in gene_list:
        sys.stdout.write(gene)
        for sample in slist:
            sys.stdout.write('\t' + str(master_dict[gene][sample]))
        print
    sys.stderr.write(date_time() + 'Fin!\n')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='SNP annotation.')
    parser.add_argument('-l', '--sample_list', action='store', dest='sample_list', help='list of sample names')
    parser.add_argument('-s', '--suffix', action='store', dest='suffix', help='Suffix of input files')
    parser.add_argument('-b', '--bed_ref', action='store', dest='bed_ref', help='Reference bed file')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_list, suffix, bed) = (inputs.sample_list, inputs.suffix, inputs.bed_ref)
    calc_coverage(sample_list, suffix, bed)
