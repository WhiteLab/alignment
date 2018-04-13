#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
import re
import os
from utility.job_manager import job_manager
import json
from math import log


def get_genes(bed):
    temp = {}
    gdict = {}
    for line in open(bed):
        info = line.rstrip('\n').split('\t')
        m = re.match('(\S+)_chr.*', info[-1])
        if m.group(1) not in temp:
            temp[m.group(1)] = 1
            gdict[m.group(1)] = info[0]
    return gdict


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['refs']['project_dir'], config_data['refs']['project'], config_data['tools']['bedtools'], \
           config_data['refs']['analysis'], config_data['refs']['capture'], config_data['params']['threads']


def get_bam_name(bnid, project_dir, project, align_dir):
    bam_dir = project_dir + project + '/' + align_dir + '/' + bnid + '/BAM/'
    bam = bam_dir + bnid + '.merged.final.bam'
    bai = bam_dir + bnid + '.merged.final.bai'
    f = 0
    if not os.path.isfile(bam):
        sys.stderr.write(date_time() + 'Bam not found in ' + bam_dir + '\n')
        f = 1
        return f, bam, bai
    if not os.path.isfile(bai):
        bai = bam_dir + bnid + '.merged.final.bam.bai'
        if not os.path.isfile(bai):
            sys.stderr.write(date_time() + 'Bam index file for ' + bnid + ' not found!  Please index first\n')
            f = 1
            return f, bam, bai

    return f, bam, bai


def get_gene_counts(ct_dict, in_dir, tier, bnid, suffix):
    for entry in open(in_dir + bnid + suffix):
        data = entry.rstrip('\n').split('\t')
        g = re.match(r'(\S+)_chr', data[3])
        ct_dict[bnid][tier]['TOTAL'] += float(data[4])
        ct_dict[bnid][tier][g.group(1)] += int(data[4])


def calc_tn_cov_ratios(cnv_dir, tum, norm, t1_genes, t2_genes, t1_suffix, t2_suffix):
    pair = tum + '_' + norm
    sys.stderr.write(date_time() + 'Collapsing coverage for ' + pair + '\n')
    out = open(cnv_dir + pair + '_cnv_estimate.txt', 'w')
    out.write('CHROM\tGENE\tTier\tTum Read ct\tNorm Read ct\tT/N ratio\tlog2 ratio\n')

    cur = {tum: {}, norm: {}}
    cur[tum]['t1'] = {key: 0 for key in t1_genes.keys()}
    cur[tum]['t1']['TOTAL'] = 0
    cur[tum]['t2'] = {key: 0 for key in t2_genes.keys()}
    cur[tum]['t2']['TOTAL'] = 0
    cur[norm]['t1'] = {key: 0 for key in t1_genes.keys()}
    cur[norm]['t1']['TOTAL'] = 0
    cur[norm]['t2'] = {key: 0 for key in t2_genes.keys()}
    cur[norm]['t2']['TOTAL'] = 0
    get_gene_counts(cur, cnv_dir, 't1', tum, t1_suffix)
    get_gene_counts(cur, cnv_dir, 't2', tum, t2_suffix)
    get_gene_counts(cur, cnv_dir, 't1', norm, t1_suffix)
    get_gene_counts(cur, cnv_dir, 't2', norm, t2_suffix)
    sys.stderr.write(date_time() + 'Calculating tier 1 gene coverage ratios for ' + pair + '\n')
    for gene in t1_genes:
        (tum_ct, norm_ct, tum_total, norm_total) = (cur[tum]['t1'][gene], cur[norm]['t1'][gene],
                                                    cur[tum]['t1']['TOTAL'], cur[norm]['t1']['TOTAL'])
        tum_rf, norm_rf, tn_ratio, log2_ratio = 0, 0, 0, 0
        if tum_total != 0:
            tum_rf = tum_ct/tum_total
        if norm_total != 0:
            norm_rf = norm_ct / norm_total
        if tum_rf > 0 and norm_rf > 0:
            tn_ratio = tum_rf / norm_rf
            log2_ratio = log(tn_ratio, 2)
        elif tum_rf == norm_rf:
            tn_ratio = 1
        elif tum_rf > 0 and norm_rf == 0:
            tn_ratio = float('inf')
            log2_ratio = float('inf')
        elif tum_rf == 0 and norm_rf > 0:
            tn_ratio = float('-inf')
            log2_ratio = float('-inf')

        out.write('\t'.join((t1_genes[gene], gene, '1', str(tum_ct), str(norm_ct), str(tn_ratio), str(log2_ratio)))
                  + '\n')
    sys.stderr.write(date_time() + 'Calculating tier 2 gene coverage ratios for ' + pair + '\n')
    for gene in t2_genes:
        (tum_ct, norm_ct, tum_total, norm_total) = (cur[tum]['t2'][gene], cur[norm]['t2'][gene],
                                                    cur[tum]['t2']['TOTAL'], cur[norm]['t2']['TOTAL'])
        tum_rf, norm_rf, tn_ratio, log2_ratio = 0, 0, 0, 0
        if tum_total != 0:
            tum_rf = tum_ct/tum_total
        if norm_total != 0:
            norm_rf = norm_ct / norm_total
        if tum_rf > 0 and norm_rf > 0:
            tn_ratio = tum_rf / norm_rf
            log2_ratio = log(tn_ratio, 2)
        elif tum_rf == norm_rf:
            tn_ratio = 1
        elif tum_rf > 0 and norm_rf == 0:
            tn_ratio = float('inf')
            log2_ratio = float('inf')
        elif tum_rf == 0 and norm_rf > 0:
            tn_ratio = float('-inf')
            log2_ratio = float('-inf')

        out.write('\t'.join((t2_genes[gene], gene, '2', str(tum_ct), str(norm_ct), str(tn_ratio), str(log2_ratio)))
                  + '\n')
    sys.stderr.write(date_time() + 'Completed estimating cnv for ' + pair + '\n')
    out.close()
    return 0


def cnv_pipe(config_file, tum_bam, norm_bam, o_flag, project2):
    (project_dir, project, bedtools, ana, bed, cores) = parse_config(config_file)
    job_list = []
    tum_id = re.match('(\d+-\d+)\.', os.path.basename(tum_bam))
    tum_id = tum_id.group(1)
    norm_id = re.match('(\d+-\d+)\.', os.path.basename(norm_bam))
    norm_id = norm_id.group(1)
    pair = tum_id + '_' + norm_id
    bed_t1 = bed.replace('.bed', '_t1.bed')
    bed_t2 = bed.replace('.bed', '_t2.bed')
    t1_genes = get_genes(bed_t1)
    t2_genes = get_genes(bed_t2)
    t1_suffix = '.t1.bedtools.coverage.txt'
    t2_suffix = '.t2.bedtools.coverage.txt'
    cnv_dir = project_dir + project + '/' + ana + '/' + pair + '/OUTPUT/'
    if not os.path.isdir(cnv_dir):
        sys.stderr.write(date_time() + 'Output path ' + cnv_dir + ' does not exist! Trying with backup project '
                         + project2 + '\n')
        cnv_dir = project_dir + project2 + '/' + ana + '/' + pair + '/OUTPUT/'
        if not os.path.isdir(cnv_dir):
            sys.stderr.write(date_time() + 'Output path ' + cnv_dir + ' does not exist! Check config and try again!\n')
            exit(1)
    clist = (tum_bam + ' -b ' + bed_t1 + ' > ', tum_bam + ' -b ' + bed_t2 + ' > ', norm_bam + ' -b ' + bed_t1
             + ' > ', norm_bam + ' -b ' + bed_t2 + ' > ')
    flist = (cnv_dir + tum_id + t1_suffix, cnv_dir + tum_id + t2_suffix, cnv_dir + norm_id + t1_suffix, cnv_dir
             + norm_id + t2_suffix)
    if o_flag == 'y':
        sys.stderr.write(date_time() + 'Overwrite yes given, creating new coverage files\n')
        for i in range(0, len(flist)):
            job_list.append(bedtools + ' coverage -abam ' + clist[i] + flist[i])
    else:
        sys.stderr.write(date_time() + 'Overwrite no given, checking for existing coverage files first\n')
        for i in range(0, len(flist)):
            if not os.path.isfile(flist[i]):
                job_list.append(bedtools + ' coverage -abam ' + clist[i] + flist[i])

    sys.stderr.write(date_time() + 'Calculating read depth for ' + pair + '\n')
    job_manager(job_list, cores)
    # process coverage files, assess cnv
    sys.stderr.write(date_time() + 'Collapsing read counts in to tiers and gene\n')
    calc_tn_cov_ratios(cnv_dir, tum_id, norm_id, t1_genes, t2_genes, t1_suffix, t2_suffix)
    sys.stderr.write(date_time() + 'CNV analysis complete for ' + pair + '\n')
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Pipeline to estimate CNV from tumor/normal bams from custom capture')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')
    parser.add_argument('-t', '--tumor', action='store', dest='tum_bam',
                        help='Tumor bam location')
    parser.add_argument('-n', '--normal', action='store', dest='norm_bam',
                        help='Normal bam location')
    parser.add_argument('-o', '--overwrite', action='store', dest='o_flag',
                        help='overwrite flag to determine whether to wrote over existing coverage files')
    parser.add_argument('-p', '--project2', action='store', dest='project2',
                        help='backup project to try for write output')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, tum_bam, norm_bam, o_flag, project2) = (inputs.config_file, inputs.tum_bam, inputs.norm_bam, inputs.o_flag,
                                                inputs.project2)
    cnv_pipe(config_file, tum_bam, norm_bam, o_flag, project2)
