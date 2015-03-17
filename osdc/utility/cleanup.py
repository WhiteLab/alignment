#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
import subprocess
import os

def cleanup(cid,vid,bid,vip):
    cname="REFS_" + bid
    sys.stderr.write(date_time() + "Unmounting " + cid + " from vm with ID " + vid + "\n")
    # need build variables to call nova successfully
    src_cmd='. /home/ubuntu/.novarc;'
    unmount_cmd="ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ubuntu@" + vip + " \"sh -s\" < /home/ubuntu/TOOLS/Scripts/utility/unmount.sh \"" + cname + "\""
    sys.stderr.write(date_time() + unmount_cmd + "\n")
    subprocess.call(unmount_cmd,shell=True)
    detach_vm=src_cmd+"nova volume-detach " + vid + " " + cid
    sys.stderr.write(date_time() + detach_vm + "\n")
    subprocess.call(detach_vm,shell=True)
    sleep_cmd='sleep 30s'
    subprocess.call(sleep_cmd,shell=True)
    delete_vol=src_cmd+"cinder delete " + cid
    sys.stderr.write(date_time() + "Deleting cinder volume " + cname + "with id " + cid + "\n")
    subprocess.call(delete_vol,shell=True)
    delete_vm=src_cmd + "nova delete " + vid
    sys.stderr.write(date_time() + "Deleting vm with id " + vid + "\n")
    subprocess.call(delete_vm,shell=True)
if __name__ == "__main__":
    import argparse

    parser=argparse.ArgumentParser(description='Attaches cinder volume with references to existing vm')
    parser.add_argument('-cid','--cinder-id',action='store',dest='cid',help='ID of attached cinder volume.  Use cinder to find')
    parser.add_argument('-vid','--virt-mach',action='store',dest='vid',help='Virtual machine id.  Use Nova to find')
    parser.add_argument('-id','--BID',action='store',dest='bid',help='Bionimbus project id')
    parser.add_argument('-ip','--ip_add',action='store',dest='vip',help='VM IP address')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    cid=inputs.cid
    vid=inputs.vid
    bid=inputs.bid
    vip=inputs.vip
    cleanup(cid,vid,bid,vip)
