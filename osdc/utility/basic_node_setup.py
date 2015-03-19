#!/usr/bin/python
# written by Miguel Brown 2015-Feb-17. Sets up a pipeline vm, creates cinder blocks, and attaches

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
import os
import re
import json
import argparse
from date_time import date_time
import subprocess
from setup_vm import setup_vm
from attach_cinder import attach_cinder

class Setup():

    def __init__(self, bid,json_config,wait):
        self.json_config = json_config
        self.wait=wait
        self.bid=bid
        self.parse_config()
        
    def parse_config(self):
        self.config_data = json.loads(open(self.json_config, 'r').read())
        (self.vm_id,self.flavor,self.key)=(self.config_data['vm_image']['id'],self.config_data['vm_image']['flavor'],self.config_data['vm_image']['key'])
        self.create_vm()

    def create_vm(self):
        (flag,vid,vip)=setup_vm(self.bid,self.vm_id,self.flavor,self.key,self.wait)
        vm = "vm_pipe_" + self.bid
        if flag==1:
            sid=self.config_data['cinder_ref']['id']
            size=self.config_data['cinder_ref']['size']
            sys.stderr.write("Attaching cinder volume...\n")
            vol=attach_cinder(sid,vid,self.bid,size,vip,self.wait)
            sys.stderr.write("Volume attached.  Go play and have fun!\n")
        else:
            sys.stderr.write("VM setup time out for " + vm + "Check connection settings or increase wait time and try again\n")


def main():
    parser=argparse.ArgumentParser(description='VM spawner for pipeline processes.  Sets up vm for sample analysis and attaches sotrage references.  BE SURE TO DOUBLE CHECK VM AND CINDER IDS OF CONFIGURATION FILES BEFORE PROCEEDING!')
    parser.add_argument('-id','--BID',action='store',dest='bid',help='Project Bionimbus ID')
    parser.add_argument('-j','--json',action='store',dest='config_file',help='JSON config file with snapshot ids and set up params')
    parser.add_argument('-w','--wait',action='store',dest='wait',help='Wait time before giving up on spawning an image.  Recommended value 600 (in seconds)')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()

    bid=inputs.bid
    config_file=inputs.config_file
    wait=inputs.wait                    
    setup=Setup(bid,config_file,wait)
if __name__ == '__main__':
  main()
