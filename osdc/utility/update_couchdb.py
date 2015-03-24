#!/usr/bin/python
import sys
import re
from date_time import date_time
from subprocess import call
from subprocess import check_output
import subprocess

def update_couchdb(fn):
    #Get uuid command
    setproxy='export http_proxy=http://cloud-proxy:3128; export https_proxy=http://cloud-proxy:3128; export no_proxy="rados-bionimbus-pdc.opensciencedatacloud.org";'
    server='https://128.135.219.167:6984'
    user='mbrown'
    pw='Jah6Eeji'
    db='igsb_qc'
    fh=open(fn,'r')
    for obj in fh:
        obj = obj.rstrip('\n')
        get_uuid=setproxy +  'curl -X GET ' + server + '/_uuids -k;'
        sys.stderr.write(date_time() + get_uuid + '\n')
        uuid_out=check_output(get_uuid,shell=True)
        m=re.findall('\"(\w+)\"',uuid_out)
        uuid=m[1]
        # typical response: {"uuids":["24ec4b43cfe304ff4709e76f7400074d"]}
        # reset proxies
        curl='curl -X PUT -d @' + obj + ' "' + server + '/' + db + '/' + uuid + '" -H "Content-Type: application/json" -k -u "' + user + ':' + pw + '"'
        couch_cmd=setproxy + curl
        # get response
        sys.stderr.write(date_time() + couch_cmd + '\n')
        result=check_output(couch_cmd,shell=True)
        result = result.rstrip('\n')
        sys.stdout.write(obj + '\t' + result + '\n')
        if result == 1:
            sys.stderr.write(date_time() + 'Database update failed for qc stats.  Check connection')
            exit(1)
    return 0
if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Update couch db with qc stats using a json object list file')
    parser.add_argument('-f','--file',action='store',dest='fn',help='qc_stats.json document list')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    fn=(inputs.fn)
    update_couchdb(fn)
