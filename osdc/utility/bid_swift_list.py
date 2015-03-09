#!/usr/bin/python
import sys
import re
from date_time import date_time
from subprocess import call
from subprocess import check_output
import subprocess 

def bid_swift_list(obj,cont,blist):
    src_cmd=". /home/ubuntu/.novarc;"
    fh=open(blist,'r')
    for bid in fh:
        swift_cmd=src_cmd + "swift list " + obj + " --prefix " + cont + "/" + bid 
        sys.stderr.write(date_time() + swift_cmd + "\n")
        try:
            check=call(swift_cmd,shell=True)
        except:
            sys.stderr.write(date_time() + "Lising of " + bid + ' of ' + cont + " from " + obj + " failed\n")

    return 0
if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Helper script to get all files associated with a BID')
    parser.add_argument('-o','--object',action='store',dest='obj',help='Swift object name, i.e. PANCAN')
    parser.add_argument('-c','--container',action='store',dest='cont',help='Swift container prefix, i.e. RAW/2015-1234')
    parser.add_argument('-l','--list',action='store',dest='blist',help='Bionimbus ID list, one per line')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (obj,cont,blist)=(inputs.obj,inputs.cont,inputs.blist)
    bid_swift_list(obj,cont,blist)