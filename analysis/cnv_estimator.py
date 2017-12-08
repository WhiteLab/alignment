#!/usr/bin/env python

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.date_time import date_time
import re
import subprocess
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
    return config_data['refs']['cont'], config_data['refs']['obj'], config_data['tools']['bedtools'], \
           config_data['refs']['analysis'], config_data['refs']['annotation'], config_data['refs']['capture']


def get_bam_name(bnid, src_cmd, cont, obj):
    list_cmd = src_cmd + 'swift list ' + cont + ' --prefix ' + obj + '/' + bnid + '/BAM/'
    sys.stderr.write(date_time() + list_cmd + '\nGetting BAM list\n')
    flist = subprocess.check_output(list_cmd, shell=True)
    # Use to check on download status
    bam = ''
    bai = ''
    dl_cmd = ''
    for fn in re.findall('(.*)\n', flist):
        # depending on software used to index, .bai extension may follow bam
        # test = re.match('^\S+_\w*\d+\.rmdup.srt.ba[m|i]$', fn) or re.match('^\S+_\w*\d+\.rmdup.srt.bam.bai$', fn)
        test = re.search(bnid + '.merged.final.ba[m|i]$', fn) or re.search(bnid + '.merged.final.bam.bai$', fn)
        if test:
            dl_cmd += src_cmd + 'swift download ' + cont + ' ' + fn + ';'
            if fn[-3:] == 'bam':
                bam = fn
            else:
                bai = fn
    return dl_cmd, bam, bai


def get_gene_counts(ct_dict, tier, bnid, suffix):
    for entry in open(bnid + suffix):
        data = entry.rstrip('\n').split('\t')
        g = re.match(r'(\S+)_chr', data[3])
        ct_dict[bnid][tier]['TOTAL'] += float(data[4])
        ct_dict[bnid][tier][g.group(1)] += int(data[4])


def calc_tn_cov_ratios(pair_list, t1_genes, t2_genes, t1_suffix, t2_suffix):
    for pair in pair_list:
        sys.stderr.write(date_time() + 'Collapsing reads for ' + pair + '\n')
        (tum, norm) = pair.split('\t')
        pair = pair.replace('\t', '_')
        out = open(pair + '_cnv_estimate.txt', 'w')
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
        get_gene_counts(cur, 't1', tum, t1_suffix)
        get_gene_counts(cur, 't2', tum, t2_suffix)
        get_gene_counts(cur, 't1', norm, t1_suffix)
        get_gene_counts(cur, 't2', norm, t2_suffix)
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


def cnv_pipe(config_file, sample_pairs, ref_mnt, cont2):
    src_cmd = '. /home/ubuntu/.novarc;'
    (cont, obj, bedtools, ana, ann, bed) = parse_config(config_file)
    bed = ref_mnt + '/' + bed
    bed_t1 = bed.replace('.bed', '_t1.bed')
    bed_t2 = bed.replace('.bed', '_t2.bed')
    pair_list = []
    t1_genes = get_genes(bed_t1)
    t2_genes = get_genes(bed_t2)
    t1_suffix = '.t1.bedtools.coverage.txt'
    t2_suffix = '.t2.bedtools.coverage.txt'
    # calc coverage for all gene capture regions
    for pairs in open(sample_pairs):
        job_list = []
        pair_set = pairs.rstrip('\n').split('\t')
        pair_list.append('\t'.join((pair_set[1], pair_set[2])))
        (tum_dl_cmd, tum_bam, tum_bai) = get_bam_name(pair_set[1], src_cmd, cont, obj)
        sys.stderr.write(date_time() + tum_dl_cmd + '\n')
        subprocess.call(tum_dl_cmd, shell=True)
        # ensure bam is complete, otherwise try a different location
        if not (os.path.isfile(tum_bam) and os.path.getsize(tum_bam) > 0):
            sys.stderr.write(date_time() + 'Download failed for ' + pair_set[1] + ', trying backup container\n')
            (tum_dl_cmd, tum_bam, tum_bai) = get_bam_name(pair_set[1], src_cmd, cont2, obj)
            sys.stderr.write(date_time() + tum_dl_cmd + '\n')
            subprocess.call(tum_dl_cmd, shell=True)
            if not (os.path.isfile(tum_bam) and os.path.getsize(tum_bam) > 0):
                sys.stderr.write(date_time() + 'Suitable bam for ' + pair_set[1] + ' not found. Check parameters!\n')
                exit(1)
        (norm_dl_cmd, norm_bam, norm_bai) = get_bam_name(pair_set[2], src_cmd, cont, obj)
        subprocess.call(norm_dl_cmd, shell=True)
        # ensure bam is complete, otherwise try a different location
        if not (os.path.isfile(norm_bam) and os.path.getsize(norm_bam) > 0):
            sys.stderr.write(date_time() + 'Download failed for ' + pair_set[2] + ', trying backup container\n')
            (norm_dl_cmd, norm_bam, norm_bai) = get_bam_name(pair_set[2], src_cmd, cont2, obj)
            sys.stderr.write(date_time() + norm_dl_cmd + '\n')
            subprocess.call(norm_dl_cmd, shell=True)
            if not (os.path.isfile(norm_bam) and os.path.getsize(norm_bam) > 0):
                sys.stderr.write(date_time() + 'Suitable bam for ' + pair_set[2] + ' not found. Check parameters!\n')
                exit(1)

        job_list.append(bedtools + ' coverage -abam ' + tum_bam + ' -b ' + bed_t1 + ' > ' + pair_set[1]
                        + t1_suffix)
        job_list.append(bedtools + ' coverage -abam ' + tum_bam + ' -b ' + bed_t2 + ' > ' + pair_set[1]
                        + t2_suffix)
        job_list.append(bedtools + ' coverage -abam ' + norm_bam + ' -b ' + bed_t1 + ' > ' + pair_set[2]
                        + t1_suffix)
        job_list.append(bedtools + ' coverage -abam ' + norm_bam + ' -b ' + bed_t2 + ' > ' + pair_set[2]
                        + t2_suffix)
        sys.stderr.write(date_time() + 'Calculating read depth for ' + pair_set[0] + '\n')
        job_manager(job_list, '8')
        cleanup = 'rm ' + tum_bam + ' ' + tum_bai + ' ' + norm_bam + ' ' + norm_bai + ';'
        sys.stderr.write('Removing bams for ' + pair_set[0])
        subprocess.call(cleanup, shell=True)
    # process coverage files, assess cnv
    sys.stderr.write(date_time() + 'Collapsing read counts in to tiers and gene\n')
    calc_tn_cov_ratios(pair_list, t1_genes, t2_genes, t1_suffix, t2_suffix)
    sys.stderr.write(date_time() + 'Fin!\n')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Pipeline for variant calls and annotation using mutect and snpEff')
    parser.add_argument('-sp', '--sample-pairs', action='store', dest='sample_pairs',
                        help='Tumor/normal sample pair list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')
    parser.add_argument('-r', '--reference', action='store', dest='ref_mnt',
                        help='Directory references are mounted, i.e. /mnt/cinder/REFS_XXX')
    parser.add_argument('-c', '--cont2', action='store', dest='cont2',
                        help='Backup container in case the first location does not work')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_pairs, config_file, ref_mnt, cont2) = (inputs.sample_pairs, inputs.config_file, inputs.ref_mnt,
                                                   inputs.cont2)
    cnv_pipe(config_file, sample_pairs, ref_mnt, cont2)
