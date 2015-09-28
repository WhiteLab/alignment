#!/usr/bin/env python
import os
import sys
import json
import re
import pdb
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from swiftclient import client
from hashlib import md5
from subprocess import check_output
from subprocess import call
  
def check_manifest(manifest, body):
    """
    check if a body is the same object described by the manifest

    :param manifest: the raw body of the manifest from swift
    :param body: a file like object to check against the manfiest
    """
    sys.stderr.write(date_time() + 'Checking manifest\n')
    manifest = json.loads(manifest)
    for segment in manifest:
        sys.stderr.write(date_time() +  segment['name'] + '\n')
        hasher = md5(body.read(segment['bytes']))
        sys.stderr.write(date_time() + '%s ?= %s' % (hasher.hexdigest(), segment['hash'] + '\n'))
        if hasher.hexdigest() != segment['hash']:
            #            return False
            sys.stderr.write('Not the same\n')
    #return True
    sys.stderr.write('The same\n')
 
 
def main():
    try:
        _prog, filename, container, manifest, var = sys.argv
    except ValueError:
        return "usage: prog.py <filename> <container> <manifest> <openstack variable file>"
    fh=open(var,'r')
    """
    export OS_TENANT_NAME=xxx
    export OS_USERNAME=xxx
    export OS_PASSWORD=xxx
    export OS_AUTH_URL="xxx"
    """
    for line in fh:
        line=line.rstrip('\n')
        info=line.split()
        pair=info[1].split('=')
        #pdb.set_trace()
        os.environ[pair[0]]=pair[1]
    fh.close()
    #    url, token = client.get_auth(os.environ['OS_AUTH_URL'], os.environ['OS_USERNAME'], os.environ['OS_PASSWORD'])
    # using client to get token and url doesn't seem to work, doing it the stupid way
    src_cmd='. ' + var + ';'
    deproxy='unset http_proxy; unset https_proxy;'
    swift_cmd=deproxy + src_cmd + "swift stat -v " + container + " " +  manifest
    sys.stderr.write(date_time() + swift_cmd + "\n")
    stat=check_output(swift_cmd,shell=True)
    header=re.search('URL: (\S+)\s+Auth Token: (\S+)\s+',stat)
    url=header.group(1)
    token=header.group(2)
    # subtract object and manifest from url
#    m=re.match('(.*)'+container+'\/manifest',url)
    url = url.replace('/' + container + '/' + manifest,'')
    sys.stderr.write(date_time() + 'URL: ' + url + ' token: ' + token + '\n')
    headers, body = client.get_object(url, token, container, manifest, query_string='multipart-manifest=get')
    sys.stderr.write(date_time() + 'Object information recieved\n')
    with open(filename) as f:
        is_valid = check_manifest(body, f)
 
    if is_valid:
        return 0
    else:
        return 1
 
 
if __name__ == "__main__":
    sys.exit(main())
