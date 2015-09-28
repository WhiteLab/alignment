#!/usr/bin/python

import sys
import re
import subprocess
import os
import pdb

fh=open(sys.argv[1],'r')
cont='HGAC'
for obj in fh:
    obj = obj.rstrip('\n')
    fn=os.path.basename(obj)
    jfile=open(obj,'r')
    temp=open(fn,'w')
    head=next(jfile)
    temp.write(head)
    garbage=next(jfile)
    for i in xrange(0,7,1):
        cur=next(jfile)
        cur=cur.replace('\t\t','\t')
        temp.write(cur)
#    pdb.set_trace()
    cur=next(jfile)
    cur=cur.replace('\n',',\n')
    cur=cur.replace('\t\t','\t')
    temp.write(cur)
    garbage=next(jfile)
    garbage=next(jfile)
    for i in xrange(0,12,1):
        cur=next(jfile)
        cur=cur.replace('\t\t','\t')
        temp.write(cur)
    cur=next(jfile)
    cur=cur.replace('\n',',\n')
    cur=cur.replace('\t\t','\t')
    temp.write(cur)
    garbage=next(jfile)
    garbage=next(jfile)
    for i in xrange(0,9,1):
        cur=next(jfile)
        cur=cur.replace('\t\t','\t')
        temp.write(cur)
    temp.write('}')
    jfile.close()
    temp.close()
    rename="mv " + fn + " " + obj
    subprocess.call(rename,shell=True)
fh.close()
