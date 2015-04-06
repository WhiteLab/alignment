#!/usr/bin/python
import sys
import re
from date_time import date_time
from subprocess import call
from subprocess import check_output
import subprocess

def delete_from_swift_list(cont,fn,l):
    src_cmd=". /home/ubuntu/.novarc;"
    deproxy='unset http_proxy; unset https_proxy;'
    fh=open(fn,'r')
    for obj in fh:
        obj = obj.rstrip('\n')
        if re.match('\W+',obj) or obj=='\n' or obj=='':
            sys.stderr.write(date_time() + 'Object ' + obj + ' looks malformed, skipping for safety reasons!\n' )
            continue
        if l== 'y':
            swift_cmd=deproxy + src_cmd + "swift delete --leave-segments " + cont + " " + obj + " >> dl_log.txt 2>> dl_log.txt"            
        else:
            swift_cmd=deproxy + src_cmd + "swift delete " + cont + " " + obj + " >> dl_log.txt 2>> dl_log.txt"
        sys.stderr.write(date_time() + swift_cmd + "\n")
        call(swift_cmd,shell=True)
    return 0
if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Delete from swift using object list')
    parser.add_argument('-c','--container',action='store',dest='cont',help='Swift container, i.e. PANCAN')
    parser.add_argument('-f','--file',action='store',dest='fn',help='Swift object list - text document one per line')
    parser.add_argument('-l','--leave',action='store',dest='l',help='Flag to leave segments (\'y\') for large objects.  Useful for when a large one was renamed, but the manifest with original anes stays the same')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (cont,fn,l)=(inputs.cont,inputs.fn,inputs.l)
    delete_from_swift_list(cont,fn,l)
