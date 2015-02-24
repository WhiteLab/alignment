#!/usr/bin/python
# written by Miguel Brown 2015-Feb-23. Wrapper script to loop through sequencing files and use pipeline

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/modules')
import os
import re
import argparse
from date_time import date_time
import subprocess

parser=argparse.ArgumentParser(description='Pipeline wrapper script to process multiple paired end set serially')
parser.add_argument('-f','--file',action='store',dest='fn',help='File with bionimbus ID, seqtype and sample list')
#parser.add_argument('-f2','--file2',action='store',dest='end2',help='Second fastq file')
#parser.add_argument('-t','--seqtype',action='store',dest='seqtype',help='Type of sequencing peformed.  Likely choices are genome, exome, and target for capture')
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

inputs=parser.parse_args()
fh=open(inputs.fn,'r')
src_cmd='. ~/.novarc;'
obj='PANCAN'

pipe_tool='/home/ubuntu/TOOLS/Scripts/modules/pipeline.py'

for line in fh:
    line=line.rstrip('\n')
    (bid,seqtype,samp_csv)=line.split('\t')
    contain='RAW/' + bid + '/' + bid + '_'
    # get list of paired end sequencing files for each pair
    for samp in samp_csv.split(', '):
        swift_cmd=src_cmd + 'swift list ' + obj + ' --prefix ' + contain + samp
        sys.stderr.write(date_time() + 'Getting sequece files for sample ' + samp + '\n' + swift_cmd + '\n')
        contents=subprocess.check_output(swift_cmd,shell=True)
        end1=''
        end2=''
        i=0
        for seqfile in re.findall('(.*)\n',contents):
            seqfile=seqfile.rstrip('\n')
            if i==0:
                end1=seqfile
                i=i+1
            else:
                end2=seqfile
                i=0
                swift_dl_cmd=src_cmd + 'swift download ' + obj + ' ' + end1
                sys.stderr.write(date_time() + 'Getting first sequencing file\n' + swift_dl_cmd + '\n')
                subprocess.call(swift_dl_cmd,shell=True)
                swift_dl_cmd=src_cmd + 'swift download ' + obj + ' ' + end2
                sys.stderr.write(date_time() + 'Getting second sequencing file\n' + swift_dl_cmd + '\n')
                subprocess.call(swift_dl_cmd,shell=True)
                pipe_cmd=pipe_tool + ' -f1 ' + end1 + ' -f2 ' + end2 + ' -t ' + seqtype + '2>> ' 'LOGS/' + samp + '.log >> ' + samp + '.log'
                sys.stderr.write(date_time() + 'Running pipeline process ' + pipe_cmd + '\n')
                status=subprocess.call(pipe_cmd,shell=True)
                if(status != 1):
                    sys.stderr.write("Pipeline process for sample " + seqfile + " failed!\n")
                    exit(3)
        
