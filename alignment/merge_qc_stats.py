#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
import re
from date_time import date_time
from subprocess import call
from subprocess import check_output
from download_from_swift import download_from_swift
import subprocess 

def download_from_swift(cont,obj,lane_list):
    src_cmd=". /home/ubuntu/.novarc;"
    lanes=open(lane_list, 'r')
    head=''
    data=[]
    for line in lanes:
        line=line.rstrip('\n')
        (bid,seqtype,lane_csv)=line.split('\t')
        for lane in lane_csv.split(', '):
            cur=obj + '/' + bid + '/QC/' + bid + '_' + lane + '.qc_stats.txt'
            swift_cmd=src_cmd + "swift download " + cont + " --skip-identical --prefix " + cur 
            sys.stderr.write(date_time() + swift_cmd + "\n")
            try:
                check=check_output(swift_cmd,shell=True,stderr=subprocess.PIPE)
            except:
                sys.stderr.write(date_time() + "Download of " + obj + " from " + cont + " failed\n")
                exit(1)
            stat=open(cur,'r')
            head=next(stat)
            data.append(next(stat))
            stat.close()
    lanes.close()
    sys.stdout.write(head)
    for datum in data:
        sys.stdout.write(datum)
    return 0
if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Uses pipeline lane list to create a summary table of qc stats')
    parser.add_argument('-c','--container',action='store',dest='cont',help='Swift container prefix, i.e. PANCAN')
    parser.add_argument('-o','--object',action='store',dest='obj',help='Swift object name/prefix, i.e. RAW/2015-1234')
    parser.add_argument('-l','--lane_list',action='store',dest='lane_list',help='Original lane list used to run pipeline')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (cont,obj,lane_list)=(inputs.cont,inputs.obj,inputs.lane_list)
    download_from_swift(cont,obj,lane_list)
