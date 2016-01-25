#!/usr/bin/env python
'''
Uses GARVIN sample names with a project ID and a lookup table with bids to create HGAC/pipeline-friendly names and
uploads to ceph.

Usage: ./garvin2hgac.py <LUT> <fq_list> <cont> <num_threads>

<LUT>           Look up table
<fq_list>       Fastq file list
<cont>          ceph object store container name
<num_threads>   Number of cores to use to multi-thread.  Given swift already does some multithreading, max/2 probably
best

Options
-h

'''
import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from docopt import docopt

args = docopt(__doc__)

lut = {}
for line in open(args['<LUT>'], 'r'):
    ids = line.rstrip('\n').split('\t')
    lut[ids[0]] = ids[1]

job_list = []
ONE_GB = 1073741824
src_cmd = '. /home/ubuntu/.novarc;'
swift_cmd = src_cmd + 'swift upload ' + args['cont'] + ' --skip-identical -S ' + str(ONE_GB) + ' --object-name '

for fq in open(args['<fq_list>'],'r'):
    info = fq.rstrip('\n').split('_')
    proj_id = fq[2] + '_' + fq[3]
    if proj_id not in lut:
        sys.stderr.write('Could not resolve project id in file name. Skipping ' + fq)
    else:
        bid = lut[proj_id]
        end = fq[-1][1]
        new_name = 'RAW/' + bid + '/' + bid + '_' + fq[2] + '_GARVIN_0000_' + fq[0] + '_' + fq[1] + '_' + end\
                   + '_sequence.txt.gz'
        sys.stderr.write(new_name + '\tnew object name to be uploaded')
        up_cmd = src_cmd + swift_cmd + new_name + ' 2>> up.log >> up.log'
        job_list.append(up_cmd)