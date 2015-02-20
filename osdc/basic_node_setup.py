#!/usr/bin/python
# written by Miguel Brown 2015-Feb-17. Sets up a pipeline vm, creates cinder blocks, and attaches

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/modules')
import os
import re
import argparse
from date_time import date_time
import subprocess

parser=argparse.ArgumentParser(description='VM spawner for pipeline processes')
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

sys.stderr.write(date_time() + "Starting vm qc for sample set " + bid + "\n")
# need build variables to call nova successfully
src_cmd='source /home/ubuntu/.novarc'
subprocess.call(src_cmd,shell=True)
vm = "vm_pipe_" + bid
vm_boot="nova boot " + vm + " --image " + image + " --flavor 5 --key_name mb_inst"
subprocess.call(vm_boot, shell=True)

# check status of vm until finshed spawing every 30s
i=30
sleep='sleep ' + str(i) + 's'
n=i
subprocess.call(sleep, shell=True)
flag=0
test='DNAseq_spawner_dev'
vip=''
vid=''
while flag==0:
    if n > wait:
        break
    else:
        sys.stderr.write(date_time() + "Checking success of vm boot. " + str(n) + " seconds have passed\n")
        check='nova list'
        p=subprocess.check_output(check,shell=True)
        for line in re.findall('(.*)\n',p):
            line=line.rstrip('\n')
            if(re.search(vm,line)):
                line=re.sub(r"\|",r"",line)
                info=line.split()
                vid=info[0]
                vname=info[1]
                vstatus=info[2]
                sys.stderr.write('Status of ' + vname + ' is ' + vstatus + ' with id ' + vid + '\n')
                if(vstatus=="ACTIVE"):
                    flag=1
                    m=re.search("private=(.*)",info[5])
                    vip=m.group(1)
                    break
        n=(n+i)
if(flag==1):
    sys.stderr.write("VM setup for " + vm + " with IP address " + vip + " with ID " + vid + " successful\n")
else:
    sys.stderr.write("VM setup time out for " + vm + "Check connection settings or increase wait time and try again\n")
