#!/usr/bin/env python
'''
Usage: ./fastq_64_to_33_batch.py <list> <th>

Options:
-h

Arguments:
<list> fastq list
<th>   num threads

'''
import sys
import os
import subprocess
sys.path.append('/home/ubuntu/TOOLS/Scripts/alignment')
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')

from docopt import docopt
args = docopt(__doc__)

fh = open(args['<list>'])
th = args['<th>']
cmd_list = []
dir_mk = 'mkdir converted'
subprocess.call(dir_mk, shell=True)
for line in fh:
    line = line.rstrip('\n')
    cmd = '/home/ubuntu/TOOLS/Scripts/utility/fastq64_to_33.py ' + line
    cmd_list.append(cmd)
from job_manager import job_manager
job_manager(cmd_list,th)
