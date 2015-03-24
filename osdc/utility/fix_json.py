#!/usr/bin/python

import sys
import re
import subprocess
import os
import pdb

src_cmd=". /home/ubuntu/.novarc;"
#Get uuid command                                                                                                                                                                                                                                                   
setproxy='export http_proxy=http://cloud-proxy:3128; export https_proxy=http://cloud-proxy:3128; export no_proxy="rados-bionimbus-pdc.opensciencedatacloud.org";'
deproxy='unset http_proxy; unset https_proxy;'
fh=open(sys.argv[1],'r')
cont='HGAC'
for obj in fh:
    obj = obj.rstrip('\n')
    fn=os.path.basename(obj)
#    pdb.set_trace()
    get_date_cmd=deproxy + src_cmd + 'swift stat ' + cont + ' ' + obj
    data=subprocess.check_output(get_date_cmd,shell='True')
    m=re.search('Last Modified: (.*GMT)',data)
    date=m.group(1)
    # go through file, make common fixes, for special cases extra fixes, add date modified, overwrite old file
    jfile=open(obj,'r')
    temp=open(fn,'w')
    for i in xrange(0,10,1):
        line=next(jfile)
        temp.write(line)
    oops=next(jfile)
    c=re.search('-418|-422|-424',obj)
    if c:
        oops=oops.replace('\n',',\n')
    temp.write(oops)
    for i in xrange(0,10,1):
        line=next(jfile)
        temp.write(line)
    oops=next(jfile)
#    pdb.set_trace()
    oops=oops.replace("tme","me")
    temp.write(oops)
    oops=next(jfile)
    oops=oops.replace("tme","me")
    temp.write(oops)
    oops=next(jfile)
    oops=oops.replace("vat","viat")
    oops=oops.replace('\n',',\n')
    c=re.search('-418|-422|-424',obj)
    if c:
        oops=oops.replace('ion:','ion":')
    oops=oops.replace('ion":','ion": ')
    temp.write(oops)
    temp.write('\t\t"date_aligned": "' + date + '"\n')
    for line in jfile:
        temp.write(line)
    temp.close()
    jfile.close()
    rename="mv " + fn + " " + obj
    subprocess.call(rename,shell=True)
fh.close()
