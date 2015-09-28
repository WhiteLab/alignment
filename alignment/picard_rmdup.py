#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from subprocess import call
from log import log

def picard_rmdup(java_tool,picard_tool,picard_tmp,sample,log_dir,ram):
    picard_rmdup_cmd=java_tool + " -Xmx" + ram + "g -jar " + picard_tool + " MarkDuplicates CREATE_INDEX=true TMP_DIR=" + picard_tmp + " REMOVE_DUPLICATES=true ASSUME_SORTED=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500 INPUT=" + sample + ".srt.bam OUTPUT=" + sample + ".rmdup.srt.bam METRICS_FILE=" + sample + ".rmdup.srt.metrics VALIDATION_STRINGENCY=LENIENT > " + log_dir + sample + ".picard.rmdup.pe.log 2>&1"
    log(log_dir + sample + ".picard.rmdup.pe.log",date_time() + picard_rmdup_cmd + "\n")
    call(picard_rmdup_cmd,shell=True) 

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Picard tools remove duplicate module.  Removes duplicates from BAM file, run after sorting BAM.')
    parser.add_argument('-j','--java',action='store',dest='java_tool',help='Java location directory, version jdk1.7.0_45 preferred')
    parser.add_argument('-p','--picard',action='store',dest='picard_tool',help='Picard jar file location')
    parser.add_argument('-pt','--picard_temp',action='store',dest='picard_tmp',help='Picard temp folder location to create')
    parser.add_argument('-sa','--sample',action='store',dest='sample',help='Sample/project name prefix')
    parser.add_argument('-l','--log',action='store',dest='log_dir',help='LOG directory location')
    parser.add_argument('-r','--ram',action='store',dest='ram',help='RAM to use in GB.  24 suggested in standard vm')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (java_tool,picard_tool,picard_tmp,sample,log_dir,ram)=(inputs.java_tool,inputs.picard_tool,inputs.picard_tmp,inputs.sample,inputs.log_dir,inputs.ram)
    picard_rmdup(java_tool,picard_tool,picard_tmp,sample,log_dir,ram)
