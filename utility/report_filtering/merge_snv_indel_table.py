#!/usr/bin/env python

import sys
import os
import re
import pdb


def process_indel_report(pair, report, merged_tbl, banned_tup, summary, length, pos_gene, min_vaf, tn_dict, norms,
                         alt_vaf, cov):
    cur = open(report)
    next(cur)
    maf = 0.01
    # biotype = 'protein_coding'
    weak_impact = {'MODIFIER': 1, 'LOW': 1}
    (tum, norm) = pair.split('_')
    for line in cur:
        if line == '\n':
            continue
        info = line.rstrip('\n').split('\t')
        vtup = '\t'.join(info[0:4])
        chr_pos = info[0] + '_' + info[1]
        if chr_pos not in pos_gene:
            # might be repeats in reports, skip if position reported already
            pos_gene[chr_pos] = {}
            pos_gene[chr_pos]['name'] = info[6]
        if pair in pos_gene[chr_pos]:
            continue
        pos_gene[chr_pos][pair] = 1
        valid = pos_gene[chr_pos]['name'] + '-' + info[0] + '_' + info[1] + '_' + info[2] + '->' + info[3]
        vap = float(info[-1]) * 100
        if pair not in tn_dict and norm in norms and valid in norms[norm]:
            merged_tbl[valid][pair] = str(vap)
            norms[norm][valid] += 1
            continue
        if vap < alt_vaf:
            if pair not in tn_dict:
                summary = update_summary(summary, pair, 'low vaf')
                continue
            elif vap < min_vaf:
                summary = update_summary(summary, pair, 'low vaf')
                continue
        if len(info[2]) > length or len(info[3]) > length:
            summary = update_summary(summary, pair, 'indel len')
            continue
        if vtup in banned_tup:
            summary = update_summary(summary, pair, 'panel')
            continue
        if info[10] in weak_impact:
            summary = update_summary(summary, pair, 'low impact')
            continue
        if len(info[5]) > 0 and float(info[5]) > maf:
            summary = update_summary(summary, pair, 'high maf')
            continue
        if int(info[-3]) + int(info[-2]) < cov:
            #if float(info[-1]) * 100 < 20:
            summary = update_summary(summary, pair, 'low coverage')
            continue

        if valid not in merged_tbl:
            merged_tbl[valid] = {}
        if valid not in norms[norm]:
            norms[norm][valid] = 0
        merged_tbl[valid][pair] = str(vap)
        norms[norm][valid] += 1
        if pair in tn_dict:
            tn_dict[pair][valid] = 1
    cur.close()
    return merged_tbl, summary, pos_gene, tn_dict, norms


def process_snv_report(pair, report, snv_head_list, snv_head_dict, tn_ratio, pos_gene, alt_vaf, cov):

    cur = open(report)
    next(cur)
    maf = 0.01
    tn = 2
    # biotype = 'protein_coding'
    weak_impact = {'MODIFIER': 1, 'LOW': 1}
    for line in cur:
        if line == '\n':
            continue
        info = line.rstrip('\n').split('\t')
        if info[-1] == 'OFF':
            continue
        chr_pos = info[snv_head_dict['chr']] + '_' + info[snv_head_dict['pos']]
        if chr_pos not in pos_gene:
            # might be repeats in reports, skip if position reported already
            pos_gene[chr_pos] = {}
            pos_gene[chr_pos]['name'] = info[snv_head_dict['gene']]
        if pair in pos_gene[chr_pos]:
            continue
        pos_gene[chr_pos][pair] = 1
        cur_vaf = info[snv_head_dict['%_tumor_alt']].rstrip('%')

        if float(cur_vaf) < alt_vaf:
            continue
        if float(info[snv_head_dict[tn_ratio]]) <= tn:
            continue
        if info[snv_head_dict['impact']] in weak_impact:
            continue
        if len(info[snv_head_dict['ExAC_MAF']]) > 0 and float(info[snv_head_dict['ExAC_MAF']]) > maf:
            continue
        if int(info[snv_head_dict['tumor_alt_count']]) + int(info[(snv_head_dict['tumor_alt_count'] - 1)]) < cov:
            continue
        sys.stdout.write(pair + '\t' + 'SNV')
        for key in snv_head_list:
            sys.stdout.write('\t' + info[snv_head_dict[key]])
        print

    cur.close()


def filter_merge_reports(reports, panel, length, alt_vaf, cov):
    snv_suffix = '.subsitutions.vep.prioritized_impact.report.xls'
    indel_suffix = '.indels.vep.prioritized_impact.report.xls'

    head = 'Sample_pair\tTYPE\tGENE\tCHROM\tPOS\tREF\tALT\tALT_CT\tALT_PCT\tsnp ID\tExAC_MAF\tIMPACT\tEFFECT\t' \
           '\tCODON_CHANGE\tAMINO_ACID_CHANGE'
    print head
    snv_head_list = ('gene', 'chr', 'pos', 'ref', 'alt', 'tumor_alt_count', '%_tumor_alt', 'snp_ID', 'ExAC_MAF',
                     'impact', 'effect', 'codon_change', 'amino_acid_change')
    snv_head_dict = {'gene': 14, 'chr': 0, 'pos': 1, 'ref': 3, 'alt': 4, 'tumor_alt_count': 9, '%_tumor_alt': 10,
                     'snp_ID': 12, 'ExAC_MAF': 13, 'impact': 17, 'effect': 16, 'codon_change': 19,
                     'amino_acid_change': 20}
    tn_ratio = 9
    indel_head_list = ('gene', 'chr', 'pos', 'ref', 'alt', 'alt_cov', 'vaf', 'snp_ID', 'ExAC_MAF',
                     'impact', 'effect', 'codon_change', 'amino_acid_change')
    indel_head_dict = {'gene': 6, 'chr': 0, 'pos': 1, 'ref': 2, 'alt': 3, 'alt_cov': 14, 'vaf':16 , 'snp_ID': 4,
                       'ExAC_MAF': 5, 'impact': 10, 'effect': 9, 'codon_change': 11, 'amino_acid_change': 12}
    vclass = 8
    pos_gene = {}
    for report in open(reports):
        report = report.rstrip('\n')
        fn = os.path.basename(report)
        m = re.search('^(\d+-\d+)(\S+)', fn)
        (pair, suffix) = (m.group(1), m.group(2))
        if suffix == snv_suffix:
            process_snv_report(pair, report, snv_head_list, snv_head_dict, tn_ratio, pos_gene, alt_vaf, cov)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Creates merged reports and filters on vaf, impact, biotype, '
                                                 't/n ratio, panel of normals, coverage.')
    parser.add_argument('-r', '--reports', action='store', dest='reports',
                        help='List of report files.')
    parser.add_argument('-p', '--panel', action='store', dest='panel',
                        help='Panel of normals')
    parser.add_argument('-c', '--min_coverage', action='store', dest='cov',
                        help='Min number of samples to see a variant to report it')

    parser.add_argument('-l', '--length', action='store', dest='length',
                        help='Max indel length')
    parser.add_argument('-v', '--vaf', action='store', dest='vaf',
                        help='Min variant allele frequency')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (reports, panel, length, avaf, cov) = (inputs.reports, inputs.panel, inputs.length, inputs.vaf, inputs.cov)
    filter_merge_reports(reports, panel, length, avaf, cov)
