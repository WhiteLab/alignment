#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from subprocess import call
import subprocess
from log import log

def novosort_sort_pe(novosort,sample,log_dir,threads,ram):
    novosort_sort_pe_cmd='mkdir novosort_tmp;' + novosort + " --threads " + threads + " --ram " + ram  + "G --tmpdir novosort_tmp --output " + sample + ".srt.bam --index  " + sample + ".bam > " + log_dir + sample + ".novosort.sort.pe.log 2>&1"
    log(log_dir + sample + ".novosort.sort.pe.log",date_time() + novosort_sort_pe_cmd + "\n")
    f=0
    try:
        f=subprocess.call(novosort_sort_pe_cmd,shell=True)
        rm_tmp='rm -rf novosort_tmp'
        subprocess.call(rm_tmp,shell=True)
    except:
        log(log_dir + sample + ".novosort.sort.pe.log",'novosort sort failed for sample ' + sample + '\n')
        exit(1)
    return f

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='novosort tool to sort BAM module.')
    parser.add_argument('-n','--novosort',action='store',dest='novosort',help='novosort binary location')
    parser.add_argument('-sa','--sample',action='store',dest='sample',help='Sample/project name prefix')
    parser.add_argument('-l','--log',action='store',dest='log_dir',help='LOG directory location')
    parser.add_argument('-t','--threads',action='store',dest='threads',help='Number of threads to use. Assuming standard vm, 8 recommonded ')
    parser.add_argument('-r','--ram',action='store',dest='ram',help='RAM to use in GB.  24 recommended for standard vm')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (novosort,sample,log_dir,threads,ram)=(inputs.novosort,inputs.sample,inputs.log_dir,inputs.threads,inputs.ram)
    novosort_sort_pe(novosort,sample,log_dir,threads,ram)
