#!/usr/bin/python
import sys
import re
from date_time import date_time
from subprocess import call
from subprocess import check_output
import subprocess 

def upload_to_swift(cont,obj):
    ONE_GB = 1073741824
    src_cmd=". /home/ubuntu/.novarc;"
    swift_cmd=src_cmd + "swift upload " + cont + " ./ --skip-identical --object-name " + obj + " -S " + str(ONE_GB)
    sys.stderr.write(date_time() + swift_cmd + "\n")
    try:
        check=check_output(swift_cmd,shell=True,stderr=subprocess.PIPE)
    except:
        sys.stderr.write(date_time() + "Upload of " + obj + " to " + cont +  " failed\n")
        exit(1)
    return 0
if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Uploads current directory contents to specified object and container')
    parser.add_argument('-o','--object',action='store',dest='obj',help='Swift object name to upload current directory contents to.  i.e. ALIGN/2015-1234')
    parser.add_argument('-c','--container',action='store',dest='cont',help='Swfit container to upload to.  i.e. PANCAN')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (cont,obj)=(inputs.cont,inputs.obj)
    upload_to_swift(cont,obj)
