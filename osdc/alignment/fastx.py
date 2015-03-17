#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from subprocess import Popen

def fastx(fastx_tool,sample,end1,end2):

    fastx_cmd='gzip -dc ' + end1 + ' | ' + fastx_tool + ' -N -o ' + sample + '_1.qs'
    sys.stderr.write(date_time() + fastx_cmd + "\n")
    Popen(fastx_cmd,shell=True,stdin=None,stdout=None,stderr=None,close_fds=True)
    fastx_cmd='gzip -dc ' + end2 + ' | ' + fastx_tool + ' -N -o ' + sample + '_2.qs'
    sys.stderr.write(date_time() + fastx_cmd + "\n")
    Popen(fastx_cmd,shell=True,stdin=None,stdout=None,stderr=None,close_fds=True)
if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='FASTX quality stats module.  Provides quality stats for fastq file and is independent of alignment.')
    parser.add_argument('-f','--fastx',action='store',dest='fastx_tool', help='Location of fastx_quality_stats tool.  Version 0.0.13.2 preferred.')
    parser.add_argument('-sa','--sample',action='store',dest='sample',help='Sample/location name prefix')
    parser.add_argument('-f1','--file1',action='store',dest='end1',help='First of paired-end fastq file')
    parser.add_argument('-f2','--file2',action='store',dest='end2',help='Second of paired-end fastq file')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (fastx_tool,sample,end1,end2)=(inputs.fastx_tool,inputs.sample,inputs.end1,inputs.end2)
    fastx(fastx_tool,sample,end1,end2)
