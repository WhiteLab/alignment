#!/usr/bin/env python

import argparse
import os
import sys
import re
from pysam import VariantFile
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from utility.log import log
from annotation.report_tools import *


def output_highest_impact(chrom, pos, ref, alt, alt_ct, non_alt_ct, vaf, ann_list, loc_dict, out, ref_flag):
    rank = ('HIGH', 'MODERATE', 'LOW', 'MODIFIER')
    top_gene = ''
    f = 0
    # secondary hit flag to avoid excessive repeats
    f1 = 0
    # index annotations by impact rank
    rank_dict = {}
    outstring = ''
    cand_top = []
    cand_next = []
    r_flag = 99
    for ann in ann_list:
        impact = ann[loc_dict['IMPACT']]
        if impact not in rank_dict:
            rank_dict[impact] = []
        rank_dict[impact].append(ann)
    for impact in rank:
        if impact in rank_dict:
            cur_rank = rank.index(impact)
            if cur_rank < r_flag:
                r_flag = cur_rank
            for ann in rank_dict[impact]:
                # need to add coverage info for indels
                (gene, tx_id, variant_class, effect, aa_pos, aa, codon, snp_id, ExAC_MAFs, biotype) = \
                    (ann[loc_dict['SYMBOL']], ann[loc_dict['Feature']], ann[loc_dict['VARIANT_CLASS']],
                     ann[loc_dict['Consequence']], ann[loc_dict['Protein_position']], ann[loc_dict['Amino_acids']],
                     ann[loc_dict['Codons']], ann[loc_dict['Existing_variation']], ann[loc_dict['ExAC_MAF']],
                     ann[loc_dict['BIOTYPE']])
                # need to parse exac maf to get desired allele freq, not all possible
                ExAC_MAF = ''
                if len(ExAC_MAFs) > 1:
                    maf_list = ExAC_MAFs.split('&')
                    for maf in maf_list:
                        check = re.match(alt + ':(\S+)', maf)
                        if check:
                            ExAC_MAF = check.group(1)
                # Format amino acid change to be oldPOSnew
                if len(aa) > 0:
                    # if a snv or syn, just aaPOS
                    test = aa.split('/')
                    if len(test) == 1:
                        aa += str(aa_pos)
                    else:
                        aa = test[0] + str(aa_pos) + test[1]
                cur_var = '\t'.join((chrom, pos, ref, alt, snp_id, ExAC_MAF, gene, tx_id, variant_class, effect,
                                     impact, biotype, codon, aa, alt_ct, non_alt_ct, vaf)) + '\n'
                if ref_flag == 'n':
                    if f == 0:
                        top_gene = gene
                        f = 1
                        outstring += cur_var
                    if f == 1 and gene != top_gene and impact != 'MODIFIER' and f1 != 0:
                        outstring += cur_var
                        f1 = 0
                else:
                    if f < 1 and cur_rank == r_flag:
                        if gene in ref_flag and tx_id == ref_flag[gene]:
                            top_gene = gene
                            outstring += cur_var
                            f = 1
                        else:
                            top_gene = gene
                            f = 0.5
                            cand_top.append(cur_var)

                    if f > 0 and gene != top_gene and impact != 'MODIFIER' and f1 < 1 and cur_rank == (r_flag + 1):
                        if gene in ref_flag and tx_id == ref_flag[gene]:
                            outstring += cur_var
                            f1 = 1
                        else:
                            f1 = 0.5
                            cand_next.append(cur_var)
    if ref_flag != 'n':
        if f == 0.5:
            outstring += cand_top[0]
        if f1 == 0.5:
            outstring += cand_next[0]
    out.write(outstring)


def gen_report(vcf, ref_flag):
    # open out file and index counts, context, etc
    fn = os.path.basename(vcf)
    parts = fn.split('.')
    loc = 'LOGS/' + parts[0] + '.indels.vep_priority.report.log'
    log(loc, date_time() + 'Creating prioritized impact reports for ' + vcf + '\n')
    vcf_in = VariantFile(vcf)

    out = open(parts[0] + '.indels.vep.prioritized_impact.report.xls', 'w')
    desired = {'Consequence': 0, 'IMPACT': 0, 'SYMBOL': 0, 'Feature': 0, 'Protein_position': 0, 'Amino_acids': 0,
               'Codons': 0, 'Existing_variation': 0, 'ExAC_MAF': 0, 'BIOTYPE': 0, 'VARIANT_CLASS': 0}

    desc_string = vcf_in.header.info['ANN'].record['Description']
    desc_string = desc_string.lstrip('"')
    desc_string = desc_string.rstrip('"')
    desc_string = desc_string.replace('Consequence annotations from Ensembl VEP. Format: ', '')
    f_pos_list = []
    desc_list = desc_string.split('|')
    ann_size = len(desc_list)
    for i in range(0, ann_size, 1):
        if desc_list[i] in desired:
            f_pos_list.append(i)
            desired[desc_list[i]] = i
    out.write('chr\tpos\tref\talt\tsnp_ID\tExAC_MAF\tgene\ttranscript_id\tvariant_class_effect\teffect\timpact'
            '\tbiotype\tcodon_change\tamino_acid_change\talt_cov\tnon_alt_cov\tvaf\n')
    if ref_flag != 'n':
        ref_flag = create_index(ref_flag)

    for record in vcf_in.fetch():
        (chrom, pos, ref, alt, alt_ct, non_alt_ct, vaf) = (record.contig, str(record.pos), record.ref, record.alts[0],
                                str(record.info['MINCOV']), str(record.info['ALTCOV']), str(record.info['COVRATIO']))
        ann_list = [_.split('|') for _ in record.info['ANN'].split(',')]
        output_highest_impact(chrom, pos, ref, alt, alt_ct, non_alt_ct, vaf, ann_list, desired, out, ref_flag)

    out.close()
    log(loc, date_time() + 'Creating prioritized report for ' + vcf + ' complete!\n')
    return 0


def main():
    parser = argparse.ArgumentParser(
        description='parse VEP annotated output into a digestable report, prioritizing highest impact calls.')
    parser.add_argument('-v', '--vcf', action='store', dest='vcf',
                        help='VEP annotated variant file')
    parser.add_argument('-r', '--reference', action='store', dest='ref', help='Tab-separated reference table with gene '
                    'symbols and refseq + ensembl ids to standardize what transcript is used.  Use flag \'n\' to skip')

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    args = parser.parse_args()

    gen_report(args.vcf, args.ref)


if __name__ == '__main__':
    main()
