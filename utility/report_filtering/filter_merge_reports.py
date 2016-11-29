#!/usr/bin/env python
import sys
import re
import os
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
#from utility.date_time import date_time


def update_summary(summary, pair, reason):
    if pair not in summary:
        summary[pair] = {}
    if reason not in summary[pair]:
        summary[pair][reason] = 1
    else:
        summary[pair][reason] += 1
    return summary


def process_indel_report(pair, report, merged_tbl, banned_tup, summary):
    cur = open(report)
    next(cur)
    maf = 0.01
    # biotype = 'protein_coding'
    weak_impact = {'MODIFIER': 1, 'LOW': 1}
    cov = 30
    for line in cur:
        if line == '\n':
            continue
        info = line.rstrip('\n').split('\t')
        vtup = '\t'.join(info[0:4])
        if vtup in banned_tup:
            summary = update_summary(summary, pair, 'panel')
            continue
        if info[10] in weak_impact:
            summary = update_summary(summary, pair, 'low impact')
            continue
        #if info[11] != biotype:
        #    summary = update_summary(summary, pair, 'wrong bioype')
        #    continue
        if len(info[5]) > 0 and float(info[5]) > maf:
            summary = update_summary(summary, pair, 'high maf')
            continue
        if int(info[-3]) + int(info[-2]) < cov:
            summary = update_summary(summary, pair, 'low coverage')
            continue
        valid = info[6] + '-' + info[0] + '_' + info[1] + '_' + info[2] + '->' + info[3]
        if valid not in merged_tbl:
            merged_tbl[valid] = {}
        merged_tbl[valid][pair] = str(float(info[-1]) * 100)
    return merged_tbl, summary


def process_snv_report(pair, report, merged_tbl, summary):
    merged_tbl[pair] = {}
    cur = open(report)
    next(cur)
    maf = 0.01
    tn = 2
    cov = 30
    # biotype = 'protein_coding'
    weak_impact = {'MODIFIER': 1, 'LOW': 1}
    for line in cur:
        if line == '\n':
            continue
        info = line.rstrip('\n').split('\t')
        valid = info[14] + '-' + info[0] + '_' + info[1] + '_' + info[3] + '->' + info[4]
        if info[-1] == 'OFF':
            summary = update_summary(summary, pair, 'off target')
            continue
        if float(info[11]) <= tn:
            summary = update_summary(summary, pair, 'low tn ratio')
            continue
        if info[17] in weak_impact:
            summary = update_summary(summary, pair, 'low impact')
            continue
        #if info[11] != biotype:
        #    summary = update_summary(summary, pair, 'wrong bioype')
        #    continue
        if len(info[13]) > 0 and float(info[13]) > maf:
            summary = update_summary(summary, pair, 'high maf')
            continue
        if int(info[5]) + int(info[6]) < cov or int(info[8]) + int(info[9]) < cov:
            summary = update_summary(summary, pair, 'low coverage')
            continue

        vaf = info[10].rstrip('%')
        if valid not in merged_tbl:
            merged_tbl[valid] = {}
        merged_tbl[valid][pair] = vaf
    return merged_tbl, summary


def filter_merge_reports(reports, panel, num_samp):
    banned_tup = {}
    summary = {}
    pairs = []
    pair_dict = {}
    reasons = ('off target', 'low tn ratio', 'low impact', 'high maf', 'low coverage', 'panel')
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
            (merged_tbl, summary) = process_indel_report(pair, report, merged_tbl, banned_tup, summary)
        else:
            (merged_tbl, summary) = process_snv_report(pair, report, merged_tbl, summary)
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
        if len(merged_tbl[var].keys()) < int(num_samp):
            min_ct += 1
            continue
        sys.stdout.write(var)
        for pair in pairs:
            if pair in merged_tbl[var]:
                sys.stdout.write('\t' + merged_tbl[var][pair])
            else:
                sys.stdout.write('\t0')
        print
    sys.stderr.write(str(min_ct) + ' variants didn\'t meet min sample count of ' + num_samp + '\n')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Creates merged reports and filters on vaf, impact, biotype, '
                                                 't/n ratio, panel of normals, coverage.')
    parser.add_argument('-r', '--reports', action='store', dest='reports',
                        help='List of report files')
    parser.add_argument('-p', '--panel', action='store', dest='panel',
                        help='Panel of normals')
    parser.add_argument('-n', '--number_samples', action='store', dest='num_samp',
                        help='Min number of samples to see a variant to report it')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (reports, panel, num_samp) = (inputs.reports, inputs.panel, inputs.num_samp)
    filter_merge_reports(reports, panel, num_samp)
