#!/usr/bin/env python

import json
import os
import re
import subprocess
import sys
import time
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.date_time import date_time
from utility.log import log


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['params']['cc_vsn'], int(config_data['params']['genome_size']), config_data['params']['ranges']


def print_header_cc(tbl, ccvsn, ranges):
    tbl.write(
        "BionimbusID\tDate\tMachine\tRun\tBarCode\tLane\tread_length\ttotal_reads\tpost_align_reads\tpercent_aligned" +
        "\tpost_rmdup_reads\tfraction_rmduped\ttarget_size\taligned_bp_ot\tfraction_aligned_bp_ot\t" +
        "fraction_sequenced_bp_ot\tmedian_insert_size\tmedian_absolute_deviation\tmean_insert_size\t" +
        "insert_standard_deviation\tdate_aligned\t" + ccvsn + "_t1_average_x_coverage\t" + ccvsn +
        "_t1_%covered_at_average\t")
    for r in ranges:
        r *= 100
        r = str(int(r))
        tbl.write(ccvsn + "_t1_" + r + "%\t")
    tbl.write(
        ccvsn + "_t1_zero_cov_bp\t" + ccvsn + "_t1_%zero_cov\t" + ccvsn + "_t2_average_x_coverage\t" + ccvsn +
        "_t2_%covered_at_average\t")
    for r in ranges:
        r *= 100
        r = str(int(r))
        tbl.write(ccvsn + "_t2_" + r + "%\t")
    tbl.write(ccvsn + "_t2_zero_cov_bp\t" + ccvsn + "_t2_%zero_cov\n")


def parseFS(FS):
    fh = open(FS, 'r')
    line = next(fh)
    line = line.rstrip('\n')
    rd_ct = re.search('^(\d+)\s', line)
    line = next(fh)
    line = next(fh)
    line = line.rstrip('\n')
    mapped = re.search('^(\d+).*\((\d+\.\d+)', line)
    fh.close()
    return rd_ct.group(1), mapped.group(1), mapped.group(2)


def parseINS(INS):
    fh = open(INS, 'r')
    for i in xrange(0, 7, 1):
        skip = next(fh)
    line = next(fh)
    line = line.rstrip('\n')
    stats = line.split('\t')
    fh.close()
    return stats[0], stats[1], stats[4], stats[5]


def parseCoverage(CF, ranges):
    fh = open(CF, 'r')
    # if run twice will end up being tier-specific
    cvg = {}
    zero_cov = next(fh)
    zero_cov.rstrip('\n')
    zero_attr = zero_cov.split('\t')
    temp = [zero_attr]
    target = int(zero_attr[3])
    zero_bp = int(zero_attr[2])
    zero_ratio = float(zero_attr[4])
    run_ratio = 1 - float(zero_attr[4])
    run_total = 0
    aln_bp = 0
    for ratio in ranges:
        cvg[ratio] = {}
        cvg[ratio]['xcov'] = 0
        cvg[ratio]['bases'] = 0
    for line in fh:
        line = line.rstrip('\n')
        attr = line.split('\t')
        temp.append(attr)
        aln_bp += (int(attr[1]) * int(attr[2]))
        run_ratio -= float(attr[4])
        run_total += int(attr[2])
        for ratio in ranges:
            if run_ratio <= ratio and cvg[ratio]['bases'] == 0:
                cvg[ratio]['bases'] = run_total
                cvg[ratio]['xcov'] = int(attr[1])
    fh.close()
    avg_cov = aln_bp / run_total
    avg_ratio = 0.0
    for i in xrange(0, len(temp), 1):
        avg_ratio += float(temp[i][4])
        # pdb.set_trace()
        if int(temp[i][1]) >= avg_cov:
            avg_ratio = 1 - avg_ratio
            break
    return target, aln_bp, run_total, cvg, zero_bp, zero_ratio, avg_cov, avg_ratio


def parse_qc(config_file, sample, cflag):
    log_dir = './'
    if os.path.isdir('LOGS'):
        log_dir = 'LOGS/'
    loc = log_dir + sample + '.qc_stats.log'
    # set up input files and variables - .hist, .flagstats, .metrics, .qs
    insert = sample + '.insert_metrics.hist'
    rawFlag = sample + '.srt.bam.flagstats'
    if os.path.isfile(rawFlag) == False:
        rawFlag = sample + '.bam.flagstats'

    rmdupFlag = sample + '.rmdup.srt.bam.flagstats'
    qs = sample + '_1.qs'
    tbl = open(sample + '.qc_stats.txt', 'w')
    json_out = open(sample + '.qc_stats.json', 'w')
    # get read length using line count of quality stats file
    fq_rd_len = subprocess.check_output('wc -l ' + qs, shell=True)
    wc_res = fq_rd_len.split()
    rd_len = int(wc_res[0]) - 1
    date_aligned = time.strftime("%c")
    (ccvsn, genome_size, range_list) = parse_config(config_file)
    ranges = range_list.split(',')
    ranges = map(float, ranges)
    # list containing read group info
    RG = sample.split('_')
    log(loc, date_time() + 'Parsing flagstats file ' + rawFlag + '\n')
    (tot_rds, post_aligned_reads, frac_aligned) = parseFS(rawFlag)
    log(loc, date_time() + 'Parsing flagstats file ' + rmdupFlag + '\n')
    (rmduped_reads, rmduped_aligned, dummy) = parseFS(rmdupFlag)
    frac_rmduped_reads = 1 - (float(rmduped_aligned) / float(post_aligned_reads))
    # estimate total base pairs covered
    all_bp_est = int(rd_len) * int(rmduped_reads)
    (target_size, aligned_bp_ot, num_bp, fraction_aligned_bp_ot, fraction_sequenced_ot) = (0, 0, 0.0, 0, 0.0)
    log(loc, date_time() + 'Parsing read insert size file ' + insert + '\n')
    (median_insert_size, median_absolute_deviation, mean_insert_size, insert_standard_deviation) = parseINS(insert)
    # output changes based on capture or whole genome
    if cflag == 'n':
        log(loc, date_time() + 'Custom capture flag given, computing stats for tier 1 and tier 2\n')
        bed_hist1 = sample + '.capture_t1.hist'
        bed_hist2 = sample + '.capture_t2.hist'
        print_header_cc(tbl, ccvsn, ranges)
        tbl.write('\t'.join(RG) + '\t' + str(rd_len))
        # cov is a dictionary with the per -ratio stats
        log(loc, date_time() + 'Parsing bedtools coverage file ' + bed_hist1 + '\n')
        (
            t1_ts, t1_aln_bp_ot, t1_num_bp, t1_cov, t1_zero_cov, t1_pct_zero_cov, t1_avg_cov,
            t1_avg_ratio) = parseCoverage(
            bed_hist1, ranges)
        target_size += t1_ts
        aligned_bp_ot += t1_aln_bp_ot
        num_bp += t1_num_bp
        log(loc, date_time() + 'Parsing bedtools coverage file ' + bed_hist2 + '\n')
        (
            t2_ts, t2_aln_bp_ot, t2_num_bp, t2_cov, t2_zero_cov, t2_pct_zero_cov, t2_avg_cov,
            t2_avg_ratio) = parseCoverage(
            bed_hist2, ranges)
        target_size += t2_ts
        aligned_bp_ot += t2_aln_bp_ot
        num_bp += t2_num_bp
        fraction_aligned_bp_ot = float(num_bp) / float(target_size)
        fraction_sequenced_ot = float(aligned_bp_ot) / float(all_bp_est)
        # store items to print as a list so that they can be cast as a string in a smoother statement
        to_print = (
            tot_rds, post_aligned_reads, frac_aligned, rmduped_reads, frac_rmduped_reads, target_size, aligned_bp_ot,
            fraction_aligned_bp_ot, fraction_sequenced_ot, median_insert_size, median_absolute_deviation,
            mean_insert_size,
            insert_standard_deviation, date_aligned, t1_avg_cov, t1_avg_ratio)
        tbl.write('\t' + '\t'.join(str(e) for e in to_print))
        for ratio in ranges:
            tbl.write('\t' + str(t1_cov[ratio]['xcov']))
        to_print_too = (t1_zero_cov, t1_pct_zero_cov, t2_avg_cov, t2_avg_ratio)
        tbl.write('\t' + '\t'.join(str(e) for e in to_print_too))
        for ratio in ranges:
            tbl.write('\t' + str(t2_cov[ratio]['xcov']))
        tbl.write('\t' + '\t'.join((str(t2_zero_cov), str(t2_pct_zero_cov))) + '\n')
        tbl.close()
        json_dict = {'BionimbusID': RG[0], 'Date': RG[1], 'Machine': RG[2], 'Run': RG[3], 'BarCode': RG[4],
                     'Lane': RG[5], 'read_length': rd_len, 'total_reads': tot_rds,
                     'post_align_reads': post_aligned_reads, 'percent_aligned': frac_aligned,
                     'post_rmdup_reads': rmduped_reads, 'fraction_rmduped': frac_rmduped_reads,
                     'target_size': target_size, 'aligned_bp_ot': aligned_bp_ot,
                     'fraction_aligned_bp_ot': fraction_aligned_bp_ot,
                     'fraction_sequenced_bp_ot': fraction_sequenced_ot, 'median_insert_size': median_insert_size,
                     'median_absolute_deviation': median_absolute_deviation, 'mean_insert_size': mean_insert_size,
                     'insert_standard_deviation': insert_standard_deviation, 'date_aligned': date_aligned, 'coverage': {
                (ccvsn + '_capture_t1'): {'average': t1_avg_cov, '%covered_at_average': t1_avg_ratio,
                                          '90%': t1_cov[0.9]['xcov'], '50%': t1_cov[0.5]['xcov'],
                                          '10%': t1_cov[0.1]['xcov'], 'zero_cov_bp': t1_zero_cov,
                                          '%zero_cov': t1_pct_zero_cov},
                (ccvsn + '_capture_t2'): {'average': t2_avg_cov, '%covered_at_average': t2_avg_ratio,
                                          '90%': t2_cov[0.9]['xcov'], '50%': t2_cov[0.5]['xcov'],
                                          '10%': t2_cov[0.1]['xcov'], 'zero_cov_bp': t2_zero_cov,
                                          '%zero_cov': t2_pct_zero_cov}}}
        json_out.write(json.dumps(json_dict, indent=4, sort_keys=True))
        json_out.close()
    else:
        log(loc, date_time() + 'Whole genome flag given\n')
        bed_hist = sample + '.genome.hist'
        # cov is a dictionary with the per -ratio stats
        tbl.write(
            "BionimbusID\tDate\tMachine\tRun\tBarCode\tLane\tread_length\ttotal_reads\tpost_align_reads\t" +
            "percent_aligned\tpost_rmdup_reads\tfraction_rmduped\ttarget_size\taligned_bp_ot\tfraction_aligned_bp_ot" +
            "\tfraction_sequenced_bp_ot\tmedian_insert_size\tmedian_absolute_deviation\tmean_insert_size\t" +
            "insert_standard_devation\tdate_aligned\taverage_x_coverage\t%covered_at_average\t90%\t50%\t10%\t" +
            "zero_cov_bp\t%zero_cov\n")
        tbl.write('\t'.join(RG) + '\t' + str(rd_len))
        log(loc, 'Parsing bedtools coverage file ' + bed_hist + '\n')
        (target_size, aligned_bp_ot, num_bp, cov, zero_cov, pct_zero_cov, avg_cov, avg_ratio) = parseCoverage(bed_hist,
                                                                                                              ranges)
        fraction_aligned_bp_ot = float(num_bp) / float(target_size)
        fraction_sequenced_bp_ot = float(aligned_bp_ot) / float(all_bp_est)
        to_print = (
            tot_rds, post_aligned_reads, frac_aligned, rmduped_reads, frac_rmduped_reads, target_size, aligned_bp_ot,
            fraction_aligned_bp_ot, fraction_sequenced_bp_ot, median_insert_size, median_absolute_deviation,
            mean_insert_size,
            insert_standard_deviation, date_aligned, avg_cov, avg_ratio)
        tbl.write('\t' + '\t'.join(str(e) for e in to_print))
        for ratio in ranges:
            tbl.write('\t' + str(cov[ratio]['xcov']))
        tbl.write('\t' + "\t".join((str(zero_cov), str(pct_zero_cov))) + '\n')
        tbl.close()
        json_dict = {'BionimbusID': RG[0], 'Date': RG[1], 'Machine': RG[2], 'Run': RG[3], 'BarCode': RG[4],
                     'Lane': RG[5], 'read_length': rd_len, 'total_reads': tot_rds,
                     'post_align_reads': post_aligned_reads, 'percent_aligned': frac_aligned,
                     'post_rmdup_reads': rmduped_reads, 'fraction_rmduped': frac_rmduped_reads,
                     'target_size': target_size, 'aligned_bp_ot': aligned_bp_ot,
                     'fraction_aligned_bp_ot': fraction_aligned_bp_ot,
                     'fraction_sequenced_bp_ot': fraction_sequenced_bp_ot, 'median_insert_size': median_insert_size,
                     'median_absolute_deviation': median_absolute_deviation, 'mean_insert_size': mean_insert_size,
                     'insert_standard_deviation': insert_standard_deviation, 'date_aligned': date_aligned, 'coverage': {
                ccvsn: {'average': avg_cov, 'covered_at_average': avg_ratio, '90%': cov[0.9]['xcov'],
                        '50%': cov[0.5]['xcov'], '10%': cov[0.1]['xcov'], 'zero_cov_bp': zero_cov,
                        '%zero_cov': pct_zero_cov}}}
        json_out.write(json.dumps(json_dict, indent=4, sort_keys=True))
        json_out.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Coverage and algnment QC summary.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-f', '--flag', action='store', dest='cflag',
                        help='\'y\' if whole genome, \'n\' if custom capture to provide per-tier coverage stats')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample prefix')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample, cflag) = (inputs.config_file, inputs.sample, inputs.cflag)
    parse_qc(config_file, sample, cflag)
