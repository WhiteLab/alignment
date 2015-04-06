#!/usr/bin/python
import sys
import re
from date_time import date_time
from subprocess import call
from subprocess import check_output
import subprocess
import pdb

def download_from_swift_list(cont,fn):
    src_cmd=". /home/ubuntu/.novarc;"
    deproxy='unset http_proxy; unset https_proxy;'
    fh=open(fn,'r')
    for obj in fh:
        obj=obj.rstrip('\n')
        (old,new)=obj.split('\t')
        swift_cmd=deproxy + src_cmd + "swift stat -v " + cont + " " +  old
        sys.stderr.write(date_time() + swift_cmd + "\n")
        stat=check_output(swift_cmd,shell=True)
#        pdb.set_trace()
        header=re.search('URL: (\S+)\s+Auth Token: (\S+)\s+',stat)
        url=header.group(1)
        token=header.group(2)
#        sys.stderr.write(url + '\n' + token)
        cp_cmd='curl -i ' + url + ' -X  COPY -H "X-Auth-Token: ' + token + '" -H "Destination: ' + new + '"'
        sys.stderr.write(cp_cmd + '\n')
        call(cp_cmd,shell=True)

    return 0
if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Tool to take a list of swift objects from the same container and copy them server-side using curl commands')
    parser.add_argument('-c','--container',action='store',dest='cont',help='Swift container, i.e. PANCAN')
    parser.add_argument('-f','--file',action='store',dest='fn',help='Tab-separated renaming list old <tab> new.  New object must start with container name in the event that ones is trying to also switch containers')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (cont,fn)=(inputs.cont,inputs.fn)
    download_from_swift_list(cont,fn)
