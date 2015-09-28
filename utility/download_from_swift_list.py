#!/usr/bin/python
import sys
import re
from date_time import date_time
from subprocess import call
from subprocess import check_output
import subprocess

def download_from_swift_list(cont,fn):
    src_cmd=". /home/ubuntu/.novarc;"
    deproxy='unset http_proxy; unset https_proxy;'
    fh=open(fn,'r')
    for obj in fh:
        swift_cmd=deproxy + src_cmd + "swift download " + cont + " --skip-identical " + obj + " >> dl_log.txt"
        sys.stderr.write(date_time() + swift_cmd + "\n")
        call(swift_cmd,shell=True)
    return 0
if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Download from swift using container list')
    parser.add_argument('-c','--container',action='store',dest='cont',help='Swift container, i.e. PANCAN')
    parser.add_argument('-f','--file',action='store',dest='fn',help='Swift object list - text document one per line')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (cont,fn)=(inputs.cont,inputs.fn)
    download_from_swift_list(cont,fn)
