#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from subprocess import check_output
import subprocess


def skip_lines(fh, stop):
    for i in range(0, stop, 1):
        skip = next(fh)
    return 0

def process_line(fh, stop):
    list_list = []
    for i in range(0,stop,1):
        cur = next(fh)
        cur = cur.rstrip('\n').split()
        list_list.append(cur)
    return list_list

def download_from_swift(cont, obj, lane_list):
    src_cmd = ". /home/ubuntu/.novarc;"
    lanes = open(lane_list, 'r')
    head = ''
    data = []
    print ('BID\tread group\ttotal alignment pairs(ap)\t% unambiguous ap\t% ambiguous ap\t% total ap filtered' \
          '\t%total ap kept')
    for line in lanes:
        line = line.rstrip('\n')
        (bid, seqtype, lane_csv) = line.split('\t')
        for lane in lane_csv.split(', '):
            cur = obj + '/' + bid + '/QC/' + bid + '_' + lane + '.runlog.txt'
            swift_cmd = src_cmd + "swift download " + cont + " --skip-identical --prefix " + cur
            sys.stderr.write(date_time() + swift_cmd + "\n")
            try:
                check = check_output(swift_cmd, shell=True, stderr=subprocess.PIPE).encode()
            except:
                sys.stderr.write(date_time() + "Download of " + obj + " from " + cont + " failed\n")
                exit(1)
            stat = open(cur, 'r')
            skip_lines(stat, 4)
            temp = []
            group = process_line(stat, 2)
            # may need to adjust or switch to regex in a case % sign present
            unamb_pairs_pct = group[0][-1][:-1]
            amb_pairs_pct = group[1][-1][:-1]
            filt = str(100-float(unamb_pairs_pct)-float(amb_pairs_pct))
            kept = str(float(unamb_pairs_pct) + float(amb_pairs_pct))
            temp.extend((group[0][6], unamb_pairs_pct, amb_pairs_pct, filt, kept))

            print (bid + '\t' + lane + '\t' + '\t'.join(temp))
            stat.close()


    lanes.close()
    sys.stdout.write(head)
    for datum in data:
        sys.stdout.write(datum)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Uses pipeline lane list to create a summary table of read filter stats')
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
