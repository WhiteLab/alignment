#!/usr/bin/python
# written by Miguel Brown 2015-Feb-11. test run based on align.pl calls just to confirm successful installation of tools on pipeline vm

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/modules')
import os
import re
import argparse
from date_time import date_time
from fastx import fastx
from bwa_mem_pe import bwa_mem_pe
from picard_sort_pe import picard_sort_pe
from picard_rmdup import picard_rmdup
from flagstats import flagstats
from coverage import *
from subprocess import call

parser=argparse.ArgumentParser(description='DNA alignment paired-end QC pipeline')
parser.add_argument('-f1','--file1',action='store',dest='end1',help='First fastq file')
parser.add_argument('-f2','--file2',action='store',dest='end2',help='Second fastq file')

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

inputs=parser.parse_args()

end1=inputs.end1
end2=inputs.end2
sys.stderr.write(date_time() + "Starting alignment qc for paired end sample files " + end1 + " and " + end2 + "\n")
#inputs

SAMPLES={}

s=re.match('^(\S+)_1_sequence\.txt\.gz$',end1)

sample=s.group(1)
HGACID=sample.split("_")
SAMPLES[sample]={}
SAMPLES[sample]['f1']=end1
SAMPLES[sample]['f2']=end2
RGRP="@RG\\tID:" + sample + "\\tLB:" + HGACID[0] + "\\tSM:" + HGACID[0] + "\\tPL:illumina"

#tools and refs

log_dir='LOGS/'
mk_log_dir='mkdir ' + log_dir
sys.stderr.write(date_time() + 'Made log directory ' + log_dir + "\n")
call(mk_log_dir,shell=True)

fastx_tool='fastx_quality_stats'
bwa_tool='/home/ubuntu/TOOLS/bwa-0.7.8/bwa'
bwa_ref='/mnt/cinder/REFS/bwa-0.7.8/hg19.fa'
samtools_tool='/home/ubuntu/TOOLS/samtools-0.1.19/samtools'
samtools_ref='/mnt/cinder/REFS/samtools-0.1.19/hg19.fa'
java_tool='/home/ubuntu/TOOLS/jdk1.7.0_45/bin/java'
picard_tool='/home/ubuntu/TOOLS/picard/dist/picard.jar'
picard_tmp='picard_tmp'
bedtools2_tool='/home/ubuntu/TOOLS/bedtools2/bin/bedtools'
exome_bed_ref='/mnt/cinder/REFS/BED/refseq.Hs19.coding_regions.merged.bed'
genome_bed_ref='/mnt/cinder/REFS/BED/hg19_complete_sorted.bed'
capture_bed_ref='/mnt/cinder/REFS/BED/capture_panel_2.0.bed'
parse_qc_stats='/home/ubuntu/TOOLS/Scripts/parse_qc_stats.pl'

wait_flag=0

fastx(fastx_tool,sample,end1,end2) # will run independently of rest of output
bwa_mem_pe(bwa_tool,RGRP,bwa_ref,end1,end2,samtools_tool,samtools_ref,sample,log_dir) # rest won't run until completed
picard_sort_pe(java_tool,picard_tool,picard_tmp,sample,log_dir) # rest won't run until completed
picard_rmdup(java_tool,picard_tool,picard_tmp,sample,log_dir)  # rest won't run until emopleted
flagstats(samtools_tool,sample) # flag determines whether to run independently or hold up the rest of the pipe until completion
exome_coverage(bedtools2_tool,sample,exome_bed_ref,wait_flag) # flag determines whether to run independently or hold up the rest of the pipe until completion
genome_coverage(bedtools2_tool,sample,genome_bed_ref,wait_flag) # flag determines whether to run independently or hold up the rest of the pipe until completion
wait_flag=1
target_coverage(bedtools2_tool,sample,capture_bed_ref,wait_flag) # flag determines whether to run independently or hold up the rest of the pipe until completion

# check to see if last expected file was generated search for .capture.hist suffix
flist=os.listdir('./')
f=0
suffix='.capture.hist'
for fn in flist:
    if fn[-13:] == suffix:
        f=1
        break
if f==1:
    sys.stderr.write(date_time() + "Pipeline process completed!\n")
else:
    sys.stderr.write(date_time() + "File with suffix " + suffix + " is missing!  If intentional, ignore this message.  Otherwise, check logs for potential failures\n")
    exit(3)
