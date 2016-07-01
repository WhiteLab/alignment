#!/usr/bin/env python
from pysam import VariantFile
import sys
import pdb
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time


def gen_report(vcf):
    vcf_in = VariantFile(vcf)
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
    sys.stdout.write('CHROM\tPOS\tREF\tAllele\tTotal Allele Count\tTotal Position Coverage\tConsequence\tIMPACT\t'
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
                sys.stdout.write(common + cur + '\n')
                temp[cur] = 1
    return 0


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='parse VEP annotated output into a digestable report.')
    parser.add_argument('-i', '--infile', action='store', dest='infile',
                        help='VEP annotated variant file')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    vcf = inputs.infile
    gen_report(vcf)
