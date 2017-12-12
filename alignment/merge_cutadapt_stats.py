#!/usr/bin/env python3

import re
import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from subprocess import check_output
import subprocess


# dumb helper function to remove parans and retun only last part of line
def process_parens(cur):
    info = cur.rstrip('\n').split()
    data = info[-1].replace('(', '')
    data = data.replace(')', '')
    return data


def parseCUTADAPT(CUTADAPT):
    try:
        fh = open(CUTADAPT, 'r')
        flag = 0
        stats = []
        while flag == 0:
            cur = next(fh)
            if re.search('Total read', cur):
                # total read pairs
                stats.append(process_parens(cur))
                cur = next(fh)
                # r1a pct
                stats.append(process_parens(cur))
                cur = next(fh)
                # r2a pct
                stats.append(process_parens(cur))
                cur = next(fh)
                # too short
                stats.append(process_parens(cur))
                cur = next(fh)
                # too rp pass
                stats.append(process_parens(cur))
                next(fh)
                flag = 1
        tot_bp_line = next(fh)
        info = tot_bp_line.split()
        tot_bp = int(info[-2].replace(',', ''))
        # total bp
        stats.append(str(tot_bp))
        next(fh)
        next(fh)
        next(fh)
        # calculate trimmed base pers per read as a pct
        r1_qt_line = next(fh)
        info = r1_qt_line.split()
        r1_pct = round(float(info[-2].replace(',', ''))/tot_bp * 100, 2)
        #r1 trimmed
        stats.append(str(r1_pct) + '%')

        r2_qt_line = next(fh)
        info = r2_qt_line.split()
        r2_pct = round(float(info[-2].replace(',', ''))/tot_bp * 100, 2)
        # r2 trimmed
        stats.append(str(r2_pct) + '%')
        # total written
        tw = next(fh)
        stats.append(process_parens(tw))
        fh.close()
        return stats
    except:
        sys.stderr.write(date_time() + 'Unable to open/process file ' + CUTADAPT + '\n')
        exit(1)
    #return tot_pairs, r1a_pct, r2a_pct, short, rp_pass, tot_bp, r1_trim, r2_trim, bp_pass


def download_from_swift(cont, obj, lane_list):
    src_cmd = ". /home/ubuntu/.novarc;"
    lanes = open(lane_list, 'r')
    head = ''
    print('BID\tread group\ttotal starting read pairs(rp)\t% r1 w/ adapter\t% r2 w/ adapter\trp too short\t% rp passed'
          '\ttotal starting base pairs(bp)\tread1 bp trimmed\tread2 bp trimmed\t% bp written')
    for line in lanes:
        line = line.rstrip('\n')
        (bid, seqtype, lane_csv) = line.split('\t')
        for lane in lane_csv.split(', '):
            cur = obj + '/' + bid + '/LOGS/' + bid + '_' + lane + '.cutadapt.log'
            swift_cmd = src_cmd + "swift download " + cont + " --skip-identical " + cur
            sys.stderr.write(date_time() + swift_cmd + "\n")
            try:
                check = check_output(swift_cmd, shell=True, stderr=subprocess.PIPE).encode()
            except:
                sys.stderr.write(date_time() + "Download of " + obj + " from " + cont + " failed\n")
                exit(1)

            temp = parseCUTADAPT(cur)

            print (bid + '\t' + lane + '\t' + '\t'.join(temp))

    lanes.close()
    sys.stdout.write(head)

    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Uses pipeline lane list to create a summary table of cutadapt stats')
    parser.add_argument('-c', '--container', action='store', dest='cont', help='Swift container prefix, i.e. PANCAN')
    parser.add_argument('-o', '--object', action='store', dest='obj',
                        help='Swift object name/prefix, i.e. RAW/2015-1234')
    parser.add_argument('-l', '--lane_list', action='store', dest='lane_list',
                        help='Original lane list used to run pipeline')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (cont, obj, lane_list) = (inputs.cont, inputs.obj, inputs.lane_list)
    download_from_swift(cont, obj, lane_list)
