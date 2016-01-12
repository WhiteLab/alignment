#!/usr/bin/env python
'''
Summarizes output of picard Collect_OxoGMetrics

Usage: ./oxog_summary.py <file_list> <suffix>

Arguments:
<file_list> list of reports to collate
<suffix>    file suffix to remove for headings

Options
-h

'''

import sys
from docopt import docopt
import re

args = docopt(__doc__)

fl = open(args['<file_list>'], 'r')
suffix = args['<suffix>']
cntxt_list = []
tbl_info = {}
df = 0
sn_list = []
for fn in fl:
    fn = fn.rstrip('\n')
    sn = re.findall('(\S+)' + suffix, fn)
    sn_list.append(sn[0])
    cur = open(fn, 'r')
    tbl_info[sn[0]] = []
    for i in xrange(0,8,1):
        skip = next(cur)
    for line in cur:
        if line == '\n':
            break
        line = line.rstrip('\n')
        info = line.split('\t')
        if df == 0:
            cntxt_list.append(info[2])
        tbl_info[sn[0]].append(float(info[11]))
    cur.close()
    df = 1
    avg = sum(tbl_info[sn[0]])/float(len(tbl_info[sn[0]]))
    tbl_info[sn[0]].append(avg)
fl.close()

cntxt_list.append('average')
sys.stdout.write('Sample/context\t' + '\t'.join(sn_list) + '\n')
for i in xrange(0,len(cntxt_list),1):
    sys.stdout.write(cntxt_list[i])
    for sn in sn_list:
        sys.stdout.write('\t' + str(tbl_info[sn][i]))
    sys.stdout.write('\n')
