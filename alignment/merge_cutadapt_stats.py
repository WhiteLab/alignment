#!/usr/bin/env python
import re
import sys
import pdb
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from subprocess import check_output
import subprocess


def skip_lines(fh, stop):
    for i in xrange(0, stop, 1):
        skip = next(fh)
    return 0

def process_line(fh, stop):
    list_list = []
    for i in xrange(0,stop,1):
        cur = next(fh)
        cur = cur.rstrip('\n').split()
        list_list.append(cur)
    return list_list

def rm_parens(line):
    line = line.replace('(','').replace(')','')
    return line

def download_from_swift(cont, obj, lane_list):
    src_cmd = ". /home/ubuntu/.novarc;"
    lanes = open(lane_list, 'r')
    head = ''
    data = []
    print 'BID\tread group\tmin length enforced\tphred score scheme\tmin score enforced\ttotal starting read pairs(rp)' \
          '\trp too short\t% rp too short\trp passed\t% rp passed\ttotal starting base pairs(bp)\tread1 bp trimmed' \
          '\tread2 bp trimmed\t% bp trimmed\ttotal bp written\t% bp written'
    for line in lanes:
        line = line.rstrip('\n')
        (bid, seqtype, lane_csv) = line.split('\t')
        for lane in lane_csv.split(', '):
            cur = obj + '/' + bid + '/LOGS/' + bid + '_' + lane + '.cutadapt.log'
            swift_cmd = src_cmd + "swift download " + cont + " --skip-identical --prefix " + cur
            sys.stderr.write(date_time() + swift_cmd + "\n")
            try:
                check = check_output(swift_cmd, shell=True, stderr=subprocess.PIPE)
            except:
                sys.stderr.write(date_time() + "Download of " + obj + " from " + cont + " failed\n")
                exit(1)
            stat = open(cur, 'r')
            skip_lines(stat, 1)
            temp = []
            # pdb.set_trace()
            params = next(stat)
            m = re.search('-m (\d+) --quality-base=(\d+) -q (\d+)', params)
            # min len enforced, scoring scheme, min base quality)
            temp.extend([m.group(1), m.group(2), m.group(3)])
            skip_lines(stat, 7)
            # use generic funtion to parse groups on lines, take what's needed, next group
            # total reads pairs summary section
            group = process_line(stat, 5)
            # pdb.set_trace()
            temp.append(group[0][4])
            # skip next lines for now, may want in future if adapter trimming performed
            group[3][6] = rm_parens(group[3][6])
            temp.extend(group[3][5:7])
            group[4][5] = rm_parens(group[4][5])
            temp.extend(group[4][4:6])
            skip_lines(stat, 1)
            group = process_line(stat, 1)
            temp.append(group[0][3])
            skip_lines(stat, 2)
            group = process_line(stat, 3)
            group[0][3] = rm_parens(group[0][3])
            temp.extend((group[1][2], group[2][2], group[0][3]))
            group = process_line(stat, 1)
            group[0][5] = rm_parens(group[0][5])
            temp.extend((group[0][3], group[0][5]))
            print bid + '\t' + lane + '\t' + '\t'.join(temp)
            stat.close()


    lanes.close()
    sys.stdout.write(head)
    for datum in data:
        sys.stdout.write(datum)
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
