#!/usr/bin/python
import sys
import re
from date_time import date_time
from subprocess import call
from subprocess import check_output
import subprocess 

def download_from_swift(obj,cont):
    src_cmd=". /home/ubuntu/.novarc;"
    swift_cmd=src_cmd + "swift download " + obj + " --skip-identical --prefix " + cont 
    sys.stderr.write(date_time() + swift_cmd + "\n")
    try:
        check=check_output(swift_cmd,shell=True,stderr=subprocess.PIPE)
    except:
        sys.stderr.write(date_time() + "Download of " + cont + " from " + obj + " failed\n")
        exit(1)
    return 0
if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Simple download module to get files from swift.  Can use prefix or whole object name')
    parser.add_argument('-o','--object',action='store',dest='obj',help='Swift object name, i.e. PANCAN')
    parser.add_argument('-c','--container',action='store',dest='cont',help='Swift container prefix, i.e. RAW/2015-1234')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (obj,cont)=(inputs.obj,inputs.cont)
    download_from_swift(obj,cont)
