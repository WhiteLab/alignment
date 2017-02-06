#!/usr/bin/env python
import sys
import re
import os
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
import pdb


def update_summary(summary, pair, reason):
    if pair not in summary:
        summary[pair] = {}
    if reason not in summary[pair]:
        summary[pair][reason] = 1
    else:
        summary[pair][reason] += 1
    return summary


def process_indel_report(pair, report, merged_tbl, banned_tup, summary, length, pos_gene, min_vaf, tn_dict, norms,
                         alt_vaf):
    cur = open(report)
    next(cur)
    maf = 0.01
    # biotype = 'protein_coding'
    weak_impact = {'MODIFIER': 1, 'LOW': 1}
    cov = 0
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


def process_snv_report(pair, report, merged_tbl, summary, pos_gene, min_vaf, tn_dict, norms, alt_vaf):
    merged_tbl[pair] = {}
    cur = open(report)
    next(cur)
    maf = 0.01
    tn = 2
    cov = 0
    (tum, norm) = pair.split('_')
    # biotype = 'protein_coding'
    weak_impact = {'MODIFIER': 1, 'LOW': 1}
    for line in cur:
        if line == '\n':
            continue
        info = line.rstrip('\n').split('\t')
        if info[-1] == 'OFF':
            summary = update_summary(summary, pair, 'off target')
        chr_pos = info[0] + '_' + info[1]
        if chr_pos not in pos_gene:
            # might be repeats in reports, skip if position reported already
            pos_gene[chr_pos] = {}
            pos_gene[chr_pos]['name'] = info[14]
        if pair in pos_gene[chr_pos]:
            continue
        pos_gene[chr_pos][pair] = 1
        valid = pos_gene[chr_pos]['name'] + '-' + info[0] + '_' + info[1] + '_' + info[3] + '->' + info[4]
        cur_vaf = info[10].rstrip('%')
        if pair not in tn_dict and norm in norms and valid in norms[norm]:
            merged_tbl[valid][pair] = cur_vaf
            norms[norm][valid] += 1
            continue
        if float(cur_vaf) < alt_vaf:
            if pair not in tn_dict:
                summary = update_summary(summary, pair, 'low vaf')
                continue
            elif float(cur_vaf) < min_vaf:
                summary = update_summary(summary, pair, 'low vaf')
                continue
        if float(info[11]) <= tn:
            summary = update_summary(summary, pair, 'low tn ratio')
            continue
        if info[17] in weak_impact:
            summary = update_summary(summary, pair, 'low impact')
            continue
        if len(info[13]) > 0 and float(info[13]) > maf:
            summary = update_summary(summary, pair, 'high maf')
            continue
        if int(info[5]) + int(info[6]) < cov or int(info[8]) + int(info[9]) < cov:
            # if float(cur_vaf) < 20:
            summary = update_summary(summary, pair, 'low coverage')
            continue

        if valid not in norms[norm]:
            norms[norm][valid] = 0
        if valid not in merged_tbl:
            merged_tbl[valid] = {}
        merged_tbl[valid][pair] = cur_vaf
        norms[norm][valid] += 1
        if pair in tn_dict:
            #pdb.set_trace()
            tn_dict[pair][valid] = 1
    cur.close()
    return merged_tbl, summary, pos_gene, tn_dict, norms


def filter_merge_reports(reports, panel, num_samp, min_type, length, vaf, tn_string, alt_vaf):
    tn_dict = {}
    norms = {}
    for tn in tn_string.split(','):
        tn_dict[tn] = {}
        (tum, norm) = tn.split('_')
        norms[norm] = {}
    banned_tup = {}
    summary = {}
    pairs = []
    pair_dict = {}
    # dict to arbitrarily choose the first gene name it sees at a position to avoid duplicates as position
    pos_gene = {}
    indel_max = int(length)
    min_vaf = float(vaf)
    alt_vaf = float(alt_vaf)
    reasons = ('off target', 'low vaf', 'low tn ratio', 'indel len', 'low impact', 'high maf', 'low coverage', 'panel')
    for line in open(panel):
        info = line.rstrip('\n').split('\t')
        banned_tup['\t'.join(info[0:4])] = 0
    merged_tbl = {}
    for report in open(reports):
        report = report.rstrip('\n')
        fn = os.path.basename(report)
        basic = re.search('^(\d+-\d+_\d+-\d+)\.(\w+)\..*', fn)
        (pair, vtype) = (basic.group(1), basic.group(2))
        sys.stderr.write('Processing pair ' + pair + ' type ' + vtype + '\n')
        if pair not in pair_dict:
            pairs.append(pair)
            pair_dict[pair] = 1
        if vtype == 'indels':
            (merged_tbl, summary, pos_gene, tn_dict, norms) = process_indel_report(pair, report, merged_tbl, banned_tup,
                                                        summary, indel_max, pos_gene, min_vaf, tn_dict, norms, alt_vaf)
        else:
            (merged_tbl, summary, pos_gene, tn_dict, norms) = process_snv_report(pair, report, merged_tbl, summary,
                                                                            pos_gene, min_vaf, tn_dict, norms, alt_vaf)
    sum_tbl = open('reject_summary.txt', 'w')
    sum_tbl.write('Pair/reason\t' + '\t'.join(reasons) + '\n')
    for pair in pairs:
        sum_tbl.write(pair)
        for reason in reasons:
            if reason in summary[pair]:
                sum_tbl.write('\t' + str(summary[pair][reason]))
            else:
                sum_tbl.write('\t0')
        sum_tbl.write('\n')
    sum_tbl.close()
    sys.stdout.write('variant/sample\t' + '\t'.join(pairs) + '\n')
    min_ct = 0
    for var in merged_tbl:
        if min_type == 'all':
            if len(merged_tbl[var].keys()) < int(num_samp):
                min_ct += 1
                continue
        var_string = var
        out_flag = 0
        if min_type == 'all':
            out_flag = 1
        for pair in pairs:
            (tum, norm) = pair.split('_')
            if min_type == 'within':
                if (pair in tn_dict and var in tn_dict[pair]) or (var in norms[norm] and norms[norm][var]
                    >= int(num_samp)):
                    out_flag = 1
            if pair in merged_tbl[var]:
                var_string += ('\t' + merged_tbl[var][pair])
            else:
                var_string += '\t0'
        if out_flag == 1:
            print var_string
        else:
            min_ct += 1
    sys.stderr.write(str(min_ct) + ' variants didn\'t meet min sample count of ' + num_samp + '\n')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Creates merged reports and filters on vaf, impact, biotype, '
                                                 't/n ratio, panel of normals, coverage.')
    parser.add_argument('-r', '--reports', action='store', dest='reports',
                        help='List of report files.  TN reports must be listed first')
    parser.add_argument('-p', '--panel', action='store', dest='panel',
                        help='Panel of normals')
    parser.add_argument('-n', '--number_samples', action='store', dest='num_samp',
                        help='Min number of samples to see a variant to report it')
    parser.add_argument('-f', '--flag_type', action='store', dest='flag_type',
                        help='Filter min sample on \'all\', or \'within\' sample group')
    parser.add_argument('-l', '--length', action='store', dest='length',
                        help='Max indel length')
    parser.add_argument('-v', '--vaf', action='store', dest='vaf',
                        help='Min variant allele frequency')
    parser.add_argument('-z', '--vflag', action='store', dest='avaf',
                        help='non-tumor/normal vaf')
    parser.add_argument('-t', '--tumor_normal', action='store', dest='tum_norm',
                        help='Comma-separated tumor-normal pair list to fix')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (reports, panel, num_samp, min_type, length, vaf, tn_string, avaf) = (inputs.reports, inputs.panel,
            inputs.num_samp, inputs.flag_type, inputs.length, inputs.vaf, inputs.tum_norm, inputs.avaf)
    filter_merge_reports(reports, panel, num_samp, min_type, length, vaf, tn_string, avaf)
