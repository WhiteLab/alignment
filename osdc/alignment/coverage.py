#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from subprocess import Popen
from subprocess import call
from log import log

def exome_coverage(bedtools2_tool,sample,exome_bed_ref,wait_flag):
    exome_coverage_cmd=bedtools2_tool + " coverage -hist -abam " + sample + ".rmdup.srt.bam -b " + exome_bed_ref + " | grep all > " + sample + ".exome.hist" 
    sys.stderr.write(date_time() + exome_coverage_cmd + "\n")
    if wait_flag ==0:
        Popen(exome_coverage_cmd,shell=True,stdin=None,stdout=None,stderr=None,close_fds=True)
    else:
        call(exome_coverage_cmd,shell=True)
def genome_coverage(bedtools2_tool,sample,genome_bed_ref,wait_flag):
    genome_coverage_cmd=bedtools2_tool + " genomecov -ibam " + sample + ".rmdup.srt.bam -g " + genome_bed_ref + " | grep genome > " + sample + ".genome.hist"
    sys.stderr.write(date_time() + genome_coverage_cmd + "\n")
    if wait_flag==0:
        Popen(genome_coverage_cmd,shell=True,stdin=None,stdout=None,stderr=None,close_fds=True)
    else:
        call(genome_coverage_cmd,shell=True)
def capture_coverage(bedtools2_tool,sample,capture_bed_ref,wait_flag):
    capture_coverage_cmd=bedtools2_tool + " coverage -hist -abam " + sample + ".rmdup.srt.bam -b " + capture_bed_ref + " | grep all > " + sample + ".capture.hist"
    sys.stderr.write(date_time() + capture_coverage_cmd + "\n")
    if wait_flag==0:
        Popen(capture_coverage_cmd,shell=True,stdin=None,stdout=None,stderr=None,close_fds=True)
    else:
        call(capture_coverage_cmd,shell=True)
if __name__ == "__main__":
    import argparse
    import coverage
    parser=argparse.ArgumentParser(description='Bedtools coverage calculation module.  Typically run last in pipeline.  See coverage parameter.')
    parser.add_argument('-bt','--bedtools',action='store',dest='bedtools2_tool',help='Location of bedtools2 tool.')
    parser.add_argument('-sa','--sample',action='store',dest='sample',help='Sample/project name prefix')
    parser.add_argument('-c','--coverage',action='store',dest='coverage',help='Name of submodule to run.  Choose from genome, exome, capture or all.')
    parser.add_argument('-bf','--bed_file',action='store',dest='bed_file',help='Bedfile list. If running all, list as string in order format \'exome,genome,capture\'.  Else, just list the one bed file')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    bedtools2_tool=inputs.bedtools2_tool
    sample=inputs.sample
    action=inputs.coverage
    bedfile=inputs.bed_file

    if action == 'all':
        blist=bedfile.split(',')
        exome_coverage(bedtools2_tool,sample,blist[0],0)
        genome_coverage(bedtools2_tool,sample,blist[1],0)
        capture_coverage(bedtools2_tool,sample,blist[2],1)
    else:
        method=getattr(coverage,(action+'_coverage'))
        method(bedtools2_tool,sample,bedfile,1)
