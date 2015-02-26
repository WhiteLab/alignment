import sys
import re
from date_time import date_time
from subprocess import call
from subprocess import check_output
import subprocess 

def upload_to_swift(bid,sample):
    ONE_GB = 1073741824
    src_cmd=". /home/ubuntu/.novarc;"
    # get the head of the sf and move so that it can be used at the end for qc stats
    mv_sf="zcat " + sample + "_1_sequence.txt.gz | head | gzip -c > ../" + sample + "_1_sequence.txt.gz "
    call(mv_sf,shell=True)
    mv_sf="zcat " + sample + "_2_sequence.txt.gz | head | gzip -c > ../" + sample + "_2_sequence.txt.gz "
    call(mv_sf,shell=True)
    p_tmp_rm="rm -rf picard_tmp"
    call(p_tmp_rm,shell=True)
    swift_cmd=src_cmd + "swift upload PANCAN ./ --skip-identical --object-name ALIGN_TEST/" + bid + " -S " + str(ONE_GB)
    sys.stderr.write(date_time() + swift_cmd + "\n")
    try:
        check=check_output(swift_cmd,shell=True,stderr=subprocess.PIPE)
    except:
        sys.stderr.write(date_time() + "Upload for " + sample + " failed\n")
        exit(1)
    sf_rm="rm " + sample + "_1_sequence.txt.gz " + sample + "_2_sequence.txt.gz"
    call(sf_rm, shell=True)
    return 0
if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Last poart of pipeline.  Uploads results to swift and clears current volume')
    parser.add_argument('-id','--bid',action='store',dest='bid',help='Bionimbus id')
    parser.add_argument('-sa','--sample',action='store',dest='sample',help='Sample/project name prefix')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (bid,sample)=(inputs.bid,inputs.sample)    
    upload_to_swift(bid,sample)
