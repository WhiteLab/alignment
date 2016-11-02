#!/usr/bin/env python

import re
from pysam import VariantFile
import sys


def output_highest_impact(chrom, pos, ref, alt, alt_ct, tot_ct, ann_list, loc_dict, out):
    rank = ('HIGH', 'MODERATE', 'LOW', 'MODIFIER')
    top_gene = ''
    f = 0
    # index annotations by impact rank
    rank_dict = {}
    outstring = ''

    for ann in ann_list:
        impact = ann[loc_dict['IMPACT']]
        if impact not in rank_dict:
            rank_dict[impact] = []
        rank_dict[impact].append(ann)
    for impact in rank:
        if impact in rank_dict:
            for ann in rank_dict[impact]:
                # need to add coverage info for indels
                (gene, variant_class, effect, aa, codon, snp_id, ExAC_MAFs, biotype, sift, clin_sig, phred) = \
                (ann[loc_dict['SYMBOL']], ann[loc_dict['VARIANT_CLASS']], ann[loc_dict['Consequence']],
                 ann[loc_dict['Amino_acids']], ann[loc_dict['Codons']], ann[loc_dict['Existing_variation']],
                 ann[loc_dict['ExAC_MAF']], ann[loc_dict['BIOTYPE']], ann[loc_dict['SIFT']], ann[loc_dict['CLIN_SIG']],
                 ann[loc_dict['CADD_PHRED']])
                # need to parse exac maf to get desired allele freq, not all possible
                ExAC_MAF = ''
                if len(ExAC_MAFs) > 1:
                    maf_list = ExAC_MAFs.split('&')
                    for maf in maf_list:
                        check = re.match(alt + ':(\S+)', maf)
                        if check:
                            ExAC_MAF = check.group(1)
                if f == 0:
                    top_gene = gene
                    f = 1
                    outstring += '\t'.join((chrom, pos, ref, alt, alt_ct, tot_ct, gene, effect, impact, biotype, codon,
                                            aa, snp_id, variant_class, sift, ExAC_MAF, clin_sig, phred)) + '\n'
                    out.write(outstring)
                if f == 1 and gene != top_gene and impact != 'MODIFIER':
                    outstring += '\t'.join((chrom, pos, ref, alt, alt_ct, tot_ct, gene, effect, impact, biotype, codon,
                                            aa, snp_id, variant_class, sift, ExAC_MAF, clin_sig, phred)) + '\n'
                    out.write(outstring)


def gen_report(vcf, sample):
    vcf_in = VariantFile(vcf)
    out = open(sample + '.germline_pass.xls', 'w')
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
            desired[desc_list[i]] = i
    out.write('CHROM\tPOS\tREF\tAllele\tTotal Allele Count\tTotal Position Coverage\tGene\tEffect\tIMPACT\t'
                    'BIOTYPE\tCodons\tAmino_acids\tExisting_variation\tVARIANT_CLASS\tSIFT\tExAC_MAF\t'
                    'CLIN_SIG\tCADD_PHRED\n')
    for record in vcf_in.fetch():
        #pdb.set_trace()
        (chrom, pos, ref, alt, alt_ct, tot_ct) = (record.contig, str(record.pos), record.ref, record.alts[0],
                                                  str(record.info['TR']), str(record.info['TC']))
        ann_list = [_.split('|') for _ in record.info['ANN'].split(',')]
        output_highest_impact(chrom, pos, ref, alt, alt_ct, tot_ct, ann_list, desired, out)
    out.close()
    return 0


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='parse VEP annotated output into a digestable report.')
    parser.add_argument('-i', '--infile', action='store', dest='infile',
                        help='VEP annotated variant file')
    parser.add_argument('-s', '--sample', action='store', dest = 'sample', help='Sample name')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (vcf, sample) = (inputs.infile, inputs.sample)
    gen_report(vcf, sample)
