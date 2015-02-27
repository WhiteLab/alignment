import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/modules')
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
import os
import re
from date_time import date_time
from fastx import fastx
from bwa_mem_pe import bwa_mem_pe
from picard_sort_pe import picard_sort_pe
from picard_rmdup import picard_rmdup
from picard_insert_size import picard_insert_size
from flagstats import flagstats
import coverage
from subprocess import call
import subprocess
import json
import pdb

class Pipeline():
    
    def __init__(self,end1,end2,seqtype,json_config):
        self.json_config = json_config
        self.end1=end1
        self.end2=end2
        self.seqtype=seqtype
        self.parse_config()

    def log(self,string):
        log_file=open('LOGS/' + self.sample + '.pipe.log','a')
        log_file.write(string)
        log_file.close()

    def parse_config(self):
        self.config_data = json.loads(open(self.json_config, 'r').read())
        s=re.match('^(\S+)_1_sequence\.txt\.gz$',self.end1)
        self.sample=s.group(1)
        HGACID=self.sample.split("_")
        self.bid=HGACID[0]
        self.ref_mnt="/mnt/cinder/REFS_" + self.bid
        self.fastx_tool=self.config_data['tools']['fastx']
        self.bwa_tool=self.config_data['tools']['bwa']
        self.bwa_ref=self.ref_mnt + '/' + self.config_data['refs']['bwa']
        self.samtools_tool=self.config_data['tools']['samtools']
        self.samtools_ref=self.ref_mnt + '/' + self.config_data['refs']['samtools']
        self.java_tool=self.config_data['tools']['java']
        self.picard_tool=self.config_data['tools']['picard']
        self.picard_tmp='picard_tmp'
        self.bedtools2_tool=self.config_data['tools']['bedtools']
        self.bed_ref=self.ref_mnt + '/' + self.config_data['refs'][self.seqtype]
        self.obj=self.config_data['refs']['obj']
        self.cont=self.config_data['refs']['cont']
        self.pipeline()

    def pipeline(self):
        log_dir='LOGS/'
        if os.path.isdir(log_dir) == False:
            mk_log_dir='mkdir ' + log_dir
            call(mk_log_dir,shell=True)
            self.log(date_time() + 'Made log directory ' + log_dir + "\n")

        self.log(date_time() + "Starting alignment qc for paired end sample files " + self.end1 + " and " + self.end2 + "\n")
        #inputs
        
        SAMPLES={}
        SAMPLES[self.sample]={}
        SAMPLES[self.sample]['f1']=self.end1
        SAMPLES[self.sample]['f2']=self.end2
        RGRP="@RG\\tID:" + self.sample + "\\tLB:" + self.bid + "\\tSM:" + self.bid + "\\tPL:illumina"
        
        #tools and refs
    
        wait_flag=1
    
        bwa_mem_pe(self.bwa_tool,RGRP,self.bwa_ref,self.end1,self.end2,self.samtools_tool,self.samtools_ref,self.sample,log_dir) # rest won't run until completed
        fastx(self.fastx_tool,self.sample,self.end1,self.end2) # will run independently of rest of output
        picard_sort_pe(self.java_tool,self.picard_tool,self.picard_tmp,self.sample,log_dir) # rest won't run until completed
        picard_rmdup(self.java_tool,self.picard_tool,self.picard_tmp,self.sample,log_dir)  # rest won't run until emopleted
        flagstats(self.samtools_tool,self.sample) # flag determines whether to run independently or hold up the rest of the pipe until completion
        
        picard_insert_size(self.java_tool,self.picard_tool,self.sample,log_dir) # get insert size metrics. 
        #figure out which coverage method to call using seqtype
        method=getattr(coverage,(self.seqtype+'_coverage'))

        method(self.bedtools2_tool,self.sample,self.bed_ref,wait_flag) # run last since this step slowest of the last
        
        # check to see if last expected file was generated search for seqtype + .hist suffix
        flist=os.listdir('./')
        f=0
        suffix=self.seqtype+'.hist'
        for fn in flist:
            if fn==(self.sample +'.' + suffix):
                f=1
                break
        if f==1:
            # get the head of the sf and move so that it can be used at the end for qc stats                                                                                          
            mv_sf="gzip -dc " + self.end1 + " | head | gzip -c > ../" + self.end1 
            call(mv_sf,shell=True)
            mv_sf="gzip -dc " + self.end2 + " | head | gzip -c > ../" + self.end2
            call(mv_sf,shell=True)
            p_tmp_rm="rm -rf picard_tmp"
            call(p_tmp_rm,shell=True)
            
            from upload_to_swift import upload_to_swift
            cont=self.cont + "/" + self.bid
            check=upload_to_swift(self.obj,cont)
            if check==0:
                self.log(date_time() + "Pipeline process completed!  Cleaning up large files for next run\n")
                clean='rm *.bam ' + self.end1 + ' ' + self.end2
                subprocess.call(clean,shell=True)
                return 0
            else:
                self.log(date_time() + "All but file upload succeeded\n")
                return 1
        else:
            self.log(date_time() + "File with suffix " + suffix + " is missing!  If intentional, ignore this message.  Otherwise, check logs for potential failures\n")
            return 1

def main():
    import argparse
    parser=argparse.ArgumentParser(description='DNA alignment paired-end QC pipeline')
    parser.add_argument('-f1','--file1',action='store',dest='end1',help='First fastq file')
    parser.add_argument('-f2','--file2',action='store',dest='end2',help='Second fastq file')
    parser.add_argument('-t','--seqtype',action='store',dest='seqtype',help='Type of sequencing peformed.  Likely choices are genome, exome, and target for capture')
    parser.add_argument('-j','--json',action='store',dest='config_file',help='JSON config file containing tool and reference locations')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()

    end1=inputs.end1
    end2=inputs.end2
    seqtype=inputs.seqtype
    config_file=inputs.config_file
    Pipeline(end1,end2,seqtype,config_file)
if __name__ == "__main__":
    main()