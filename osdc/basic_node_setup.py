#!/usr/bin/python
# written by Miguel Brown 2015-Feb-17. Sets up a pipeline vm, creates cinder blocks, and attaches

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/modules')
import os
import re
import argparse
from date_time import date_time
import subprocess
from setup_vm import setup_vm
from attach_cinder import attach_cinder

parser=argparse.ArgumentParser(description='VM spawner for pipeline processes.  Sets up vm for sample analysis and attaches sotrage references')
parser.add_argument('-id','--BID',action='store',dest='bid',help='Project Bionimbus ID')
parser.add_argument('-im','--image',action='store',dest='image',help='Image id to spawn')
parser.add_argument('-w','--wait',action='store',dest='wait',help='Wait time before giving up on spawning an image.  Reommended value 300 (in seconds)')

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

inputs=parser.parse_args()

bid=inputs.bid
image=inputs.image
wait=inputs.wait

(flag,vid,vip)=setup_vm(bid,image,wait)
if flag==1:
    #default sid  change here for a differenet reference volume
    sid='18c83d7d-dd26-4fa3-8f35-ff378c580c27'
    sys.stderr.write("Attaching cinder volume...\n")
    vol=attach_cinder(sid,vid,bid,wait)
    sys.stderr.write("Volume attached.  Go play and have fun!\n")
else:
    sys.stderr.write("VM setup time out for " + vm + "Check connection settings or increase wait time and try again\n")
