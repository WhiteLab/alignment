#!/usr/bin/python
import sys
from date_time import date_time
from subprocess import call
import subprocess
from log import log

def novosort_sort_pe(novosort,sample,log_dir):
    novosort_sort_pe_cmd=novosort + " --threads 8 --ram 28G --output " + sample + ".srt.bam --index  " + sample + ".bam > " + log_dir + sample + ".novosort.sort.pe.log 2>&1"
    log(log_dir + sample + ".novosort.sort.pe.log",date_time() + novosort_sort_pe_cmd + "\n")
    try:
        subprocess.check_output(novosort_sort_pe_cmd,shell=True) 
    except:
        log(log_dir + sample + ".novosort.sort.pe.log",'novosort sort failed for sample ' + sample + '\n')
        exit(1)
    return 0

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='novosort tool to sort BAM module.')
    parser.add_argument('-n','--novosort',action='store',dest='novosort',help='novosort binary location')
    parser.add_argument('-sa','--sample',action='store',dest='sample',help='Sample/project name prefix')
    parser.add_argument('-l','--log',action='store',dest='log_dir',help='LOG directory location')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (novosort,sample,log_dir)=(inputs.novosort,inputs.sample,inputs.log_dir)    
    novosort_sort_pe(novosort,sample,log_dir)
