import sys
from date_time import date_time
from subprocess import call
from log import log

def picard_insert_size(java_tool,picard_tool,sample,log_dir):
    picard_insert_size_cmd=java_tool + " -Xmx2g -jar " + picard_tool + " CollectInsertSizeMetrics I=" + sample + ".rmdup.srt.bam H=" + sample + ".insert_metrics.pdf O=" + sample + ".insert_metrics.hist  > " + log_dir + sample + ".picard.insert_size.log 2>&1"
    log(log_dir + sample + ".picard.insert_size.log",date_time() + picard_insert_size_cmd + "\n")
    call(picard_insert_size_cmd,shell=True)

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Picard collect insert size metrics module.  Gathers insert size metrics, run after removing BAM duplicates.')
    parser.add_argument('-j','--java',action='store',dest='java_tool',help='Java location directory, version jdk1.7.0_45 preferred')
    parser.add_argument('-p','--picard',action='store',dest='picard_tool',help='Picard jar file location')
    parser.add_argument('-sa','--sample',action='store',dest='sample',help='Sample/project name prefix')
    parser.add_argument('-l','--log',action='store',dest='log_dir',help='LOG directory location')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (java_tool,picard_tool,sample,log_dir)=(inputs.java_tool,inputs.picard_tool,inputs.sample,inputs.log_dir)
    picard_insert_size(java_tool,picard_tool,sample,log_dir)
