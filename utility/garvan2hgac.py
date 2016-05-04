#!/usr/bin/env python
'''
Uses GARVAN sample names with a project ID and a lookup table with bids to create HGAC/pipeline-friendly names and
uploads to ceph.

Usage: ./garvan2hgac.py <LUT> <fq_list> <cont> <num_threads>

<LUT>           Look up table
<fq_list>       Fastq file list
<cont>          ceph object store container name
<num_threads>   Number of cores to use to multi-thread.  Given swift already does some multithreading, max/2 probably
best

Options
-h

'''
import sys
# import pdb
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from job_manager import job_manager
from docopt import docopt

args = docopt(__doc__)

lut = {}
for line in open(args['<LUT>'], 'r'):
    ids = line.rstrip('\n').split('\t')
    lut[ids[0]] = ids[1]

job_list = []
ONE_GB = 1073741824
src_cmd = '. /home/ubuntu/.novarc;'
swift_start = src_cmd + 'swift upload ' + args['<cont>'] + ' '
swift_end = ' --skip-identical -S ' + str(ONE_GB) + ' --object-name '

lut_out = {}

for fq in open(args['<fq_list>'], 'r'):
    fq = fq.rstrip('\n')
    info = fq.split('_')
    proj_id = info[2] + '_' + info[3]
    # pdb.set_trace()
    if proj_id not in lut:
        sys.stderr.write('Could not resolve project id in file name. Skipping ' + fq + '\n')
    else:
        bid = lut[proj_id]
        if bid not in lut_out:
            lut_out[bid] = {}
            lut_out[bid]['lane'] = []
            lut_out[bid]['orig'] = []
            lut_out[bid]['new'] = []
        end = info[-1][1]
        new_name = 'RAW/' + bid + '/' + bid + '_' + info[2] + '_GARVIN_0000_' + info[0] + '_' + info[1] + '_' + end\
                   + '_sequence.txt.gz'
        lut_out[bid]['lane'].append(info[1])
        lut_out[bid]['orig'].append(fq)
        lut_out[bid]['new'].append(new_name)

        sys.stderr.write(new_name + '\tnew object name to be uploaded')
        up_cmd = swift_start + fq + swift_end + new_name + ' 2>> up.log >> up.log'
        sys.stderr.write(up_cmd + '\n')
        job_list.append(up_cmd)
job_manager(job_list, args['<num_threads>'])

for bid in lut_out:
    sys.stdout.write(bid + '\t' + ', '.join(lut_out[bid]['orig']) + '\t' + ', '.join(lut_out[bid]['new']) + '\n')
