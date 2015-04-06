#!/usr/bin/python
import sys
from date_time import date_time
import subprocess
import os
import re

def setup_vm(bid,image,flavor,key,wait):
    sys.stderr.write(date_time() + "Starting vm qc for sample set " + bid + "\n")
    # need build variables to call nova successfully
    src_cmd='. /home/ubuntu/.novarc;'
    vm = "vm_pipe_" + bid
    vid=''
    vm_boot=src_cmd + "nova boot " + vm + " --image " + image + " --flavor " + flavor + " --key_name " + key
    sys.stderr.write(date_time() + "Booting up vm\n" + vm_boot + "\n")
    v=subprocess.check_output(vm_boot, shell=True)
    # get id of vm in the event another has the same display name
    for line in re.findall('(.*)\n',v):
        line=line.rstrip('\n')
        m=re.match('\|\s+(\S+)\s+\|\s+(\S+)\s+\|',line)
        if m:
            if m.group(1)=='id':
                vid=m.group(2)
                break
    # check status of vm until finshed spawing every 30s                                                                                                                         
    i=30
    sleep='sleep ' + str(i) + 's'
    n=i
    flag=0
    vip=''

    while flag==0:
        subprocess.call(sleep, shell=True)
        if n > int(wait):
            break
        else:
            sys.stderr.write(date_time() + "Checking success of vm boot. " + str(n) + " seconds have passed\n")
            check=src_cmd + 'nova list'
            p=subprocess.check_output(check,shell=True)
            for line in re.findall('(.*)\n',p):
                line=line.rstrip('\n')
                if(re.search(vid,line)):
                    line=re.sub(r"\|",r"",line)
                    info=line.split()
                    vname=info[1]
                    vstatus=info[2]
                    sys.stderr.write('Status of ' + vname + ' is ' + vstatus + ' with id ' + vid + '\n')
                    if(vstatus=="ACTIVE"):
                        flag=1
                        m=re.search("private=(.*)",info[5])
                        vip=m.group(1)
                        break
                    if(vstatus=="ERROR"):
                        flag=2
                        break
        n=n+i
    if(flag==1):
        # upload openstack variables from head vm
        delay='sleep 60s'
        sys.stderr.write(date_time() + 'Pausing 1 minute to give vm a chance to initialize\n')
        subprocess.call(delay,shell=True)
        nova_var='ssh-keyscan ' + vip + ' >> ~/.ssh/known_hosts;rsync /home/ubuntu/.novarc ubuntu@' + vip + ':/home/ubuntu'
        sys.stderr.write(date_time() + 'Copying openstack variables to vm\n' + nova_var + '\n')
        subprocess.call(nova_var,shell=True)
        sys.stderr.write(date_time() + "VM setup for " + vm + " with IP address " + vip + " with ID " + vid + " successful\n")
    else:
        sys.stderr.write(date_time() + "VM setup either timed out or produced ERROR for " + vm + ". Check connection settings, increase wait time and try again\n")
    return (flag,vid,vip)
                            
if __name__ == "__main__":
    import argparse

    parser=argparse.ArgumentParser(description='VM spawner for pipeline processes')
    parser.add_argument('-id','--BID',action='store',dest='bid',help='Project Bionimbus ID')
    parser.add_argument('-im','--image',action='store',dest='image',help='Image id to spawn')
    parser.add_argument('-w','--wait',action='store',dest='wait',help='Wait time before giving up on spawning an image.  Reommended value 300 (in seconds)')
    parser.add_argument('-f','--flavor',action='store',dest='flavor',help='Image \"flavor\" to spawn')
    parser.add_argument('-k','--key',action='store',dest='key',help='Image key-pair to use')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    bid=inputs.bid
    image=inputs.image
    flavor=inputs.flavor
    key=inputs.key
    wait=inputs.wait

    setup_vm(bid,image,flavor,key,wait)
