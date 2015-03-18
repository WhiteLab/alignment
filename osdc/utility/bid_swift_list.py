#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
import re
from date_time import date_time
from subprocess import call
from subprocess import check_output
import subprocess 

def bid_swift_list(cont,obj,blist):
    src_cmd=". /home/ubuntu/.novarc;"
    fh=open(blist,'r')
    for bid in fh:
        swift_cmd=src_cmd + "swift list " + cont + " --prefix " + obj + "/" + bid 
        sys.stderr.write(date_time() + swift_cmd + "\n")
        try:
            check=call(swift_cmd,shell=True)
        except:
            sys.stderr.write(date_time() + "Lising of " + bid + ' of ' + obj + " from " + cont + " failed\n")

    return 0
if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Helper script to get all files associated with a BID')
    parser.add_argument('-c','--container',action='store',dest='cont',help='Swift container name, i.e. PANCAN')
    parser.add_argument('-o','--object',action='store',dest='obj',help='Swift object prefix, i.e. RAW/2015-1234')
    parser.add_argument('-l','--list',action='store',dest='blist',help='Bionimbus ID list, one per line')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (cont,obj,blist)=(inputs.cont,inputs.obj,inputs.blist)
    bid_swift_list(cont,obj,blist)
