#!/usr/bin/env python3

import pdb
from pysam import VariantFile
import sys


def output_highest_impact(chrom, pos, ref, alt, alt_ct, tot_ct, variant_class, snp_id, ann_list, loc_dict, out):
    rank = ('HIGH', 'MODERATE', 'LOW', 'MODIFIER')
    top_gene = ''
    f = 0
    # secondary hit flag to avoid excessive repeats
    f1 = 0
    # index annotations by impact rank
    rank_dict = {}
    outstring = ''

    for ann in ann_list:
        impact = ann[loc_dict['Effect_Impact']]
        if impact not in rank_dict:
            rank_dict[impact] = []
        rank_dict[impact].append(ann)
    for impact in rank:
        if impact in rank_dict:
            for ann in rank_dict[impact]:
                # need to add coverage info for indels
                (gene, effect, aa, codon, biotype) = (ann[loc_dict['Gene_Name']], ann[loc_dict['Effect']],
                 ann[loc_dict['Amino_Acid_Change']], ann[loc_dict['Codon_Change']], ann[loc_dict['Transcript_BioType']])
                # Format amino acid change to be oldPOSnew

                cur_var = '\t'.join((chrom, pos, ref, alt, alt_ct, tot_ct, gene, effect, impact, biotype, codon,
                                            aa, snp_id, variant_class)) + '\n'
                if f == 0:
                    top_gene = gene
                    f = 1
                    outstring += cur_var
                if f == 1 and gene != top_gene and impact != 'MODIFIER' and f1 != 0:
                    outstring += cur_var
                    f1 = 1
    out.write(outstring)


def gen_report(vcf, sample):
    vcf_in = VariantFile(vcf)
    out = open(sample + '.germline_snpEff_impact_filtered_pass.xls', 'w')
    desired = {'Effect': 0, 'Effect_Impact': 0, 'Gene_Name': 0,
               'Amino_Acid_Change': 0, 'Codon_Change': 0, 'Transcript_BioType': 0, 'Existing_variation': 0}

    desc_string = vcf_in.header.info['EFF'].record['Description']
    desc_string = desc_string.replace('"', '')
    desc_string = desc_string.replace('\'', '')
    desc_string = desc_string.replace('Predicted effects for this variant.Format: ', '')
    desc_string = desc_string.replace(' ', '')
    desc_string = desc_string.replace('(','|')
    f_pos_list = []
    desc_list = desc_string.split('|')
    ann_size = len(desc_list)
    for i in range(0, ann_size, 1):
        if desc_list[i] in desired:
            f_pos_list.append(i)
            desired[desc_list[i]] = i
    out.write('CHROM\tPOS\tREF\tAllele\tTotal Allele Count\tTotal Position Coverage\tGene\t\tEffect\t'
              'IMPACT\tBIOTYPE\tCodons\tAmino_acids\tExisting_variation\tVARIANT_CLASS\n')
    for record in vcf_in.fetch():
        if 'PASS' in record.filter:

            (chrom, pos, ref, alt, alt_ct, tot_ct, existing_variation) = (record.contig, str(record.pos),
                record.ref, record.alts[0], str(record.info['TR']), str(record.info['TC']), record.id)
            if existing_variation is None:
                existing_variation = ''
            variant_class = 'SNV'
            if len(alt) > len(ref):
                variant_class = 'INSERTION'
            if len(ref) > len(alt):
                variant_class = 'DELETION'
            cur_info = record.info['EFF']
            cur_info = cur_info.replace('(','|')
            ann_list = [_.split('|') for _ in cur_info.split(',')]
            output_highest_impact(chrom, pos, ref, alt, alt_ct, tot_ct, variant_class, existing_variation, ann_list,
                                  desired, out)
    out.close()
    return 0


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='parse snpEff annotated output into a digestable report.')
    parser.add_argument('-i', '--infile', action='store', dest='infile',
                        help='snpEff annotated variant file')
    parser.add_argument('-s', '--sample', action='store', dest = 'sample', help='Sample name')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (vcf, sample) = (inputs.infile, inputs.sample)
    gen_report(vcf, sample)
