#!/usr/bin/env python

import argparse
import os
import sys
from pysam import VariantFile
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.date_time import date_time


def create_ind(out):
    res_dict = {}
    mut_out = open(out)
    next(mut_out)
    head = next(mut_out)
    head = head.rstrip('\n').split('\t')
    h_list = (2, 25, 26, 37, 38)
    for line in mut_out:
        info = line.rstrip('\n').split('\t')
        var_tup = info[0] + '\t' + info[1] + info[2] + info[3]
        res_dict[var_tup] = {}
        for i in h_list:
            res_dict[var_tup][head[i]] = info[i]
    mut_out.close()
    return res_dict


def gen_report(vcf, out, c):
    # open out file and index counts, context, etc
    (mut_dict) = create_ind(out)
    vcf_in = VariantFile(vcf)
    fn = os.path.basename(vcf)
    parts = fn.split('.')
    out = open(parts[0] + '.germline_pass.xls', 'w')
    desired = {'Consequence': 0, 'IMPACT': 0, 'SYMBOL': 0, 'Amino_acids': 0, 'Codons': 0, 'BIOTYPE': 0, 'SIFT': 0,
               'Existing_variation': 0, 'VARIANT_CLASS': 0, 'ExAC_MAF': 0, 'CLIN_SIG': 0, 'CADD_PHRED': 0}

    desc_string = vcf_in.header.info['ANN'].record['Description']
    desc_string = desc_string.lstrip('"')
    desc_string = desc_string.rstrip('"')
    desc_string = desc_string.replace('Consequence annotations from Ensembl VEP. Format: ', '')
    f_pos_list = []
    desc_list = desc_string.split('|')
    ann_size = len(desc_list)
    for i in xrange(0, ann_size, 1):
        if desc_list[i] in desired:
            f_pos_list.append(i)
            desired[desc_list[i]] = 1
    out.write('CHROM\tPOS\tREF\tAllele\tTotal Allele Count\tTotal Position Coverage\tConsequence\tIMPACT\t'
                    'SYMBOL\tBIOTYPE\tAmino_acids\tCodons\tExisting_variation\tVARIANT_CLASS\tSIFT\tExAC_MAF\t'
                    'CLIN_SIG\tCADD_PHRED\n')
    for record in vcf_in.fetch():
        #pdb.set_trace()
        common = '\t'.join((record.contig, str(record.pos), record.ref, str(record.alts), str(record.info['TR']),
                            str(record.info['TC'])))
        ann_list = [_.split('|') for _ in record.info['ANN'].split(',')]
        temp = {}
        for ann in ann_list:
            cur = ''
            for i in f_pos_list:
                cur += '\t' + ann[i]
            if cur not in temp:
                out.write(common + cur + '\n')
                temp[cur] = 1
    out.close()
    return 0




def main():
    parser = argparse.ArgumentParser(
        description='parse VEP annotated output into a digestable report, prioritizing highest impact calls.')
    parser.add_argument('-v', '--vcf', action='store', dest='vcf',
                        help='VEP annotated variant file')
    parser.add_argument('-o', '--out', action='store', dest='out',
                        help='MuTect output table - has counts for tumor/normal')
    parser.add_argument('-c', '--custom', action='store', dest='c',
                        help='bed file to mark whether hit was on or off-target. if not desired, enter \'n\' ')

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    args = parser.parse_args()

    gen_report(args.vcf, args.out, args.c)


if __name__ == '__main__':
    main()