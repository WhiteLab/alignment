import sys
from date_time import date_time
from subprocess import call

def upload_to_swift(bid,sample):
    ONE_GB = 1073741824
    setup_cmd="source /home/ubuntu/.novarc"
    call(setup_cmd, shell=True)
    sf_rm="rm " + sample + "_1_sequence.txt.gz " + sample + "_2_sequence.txt.gz"
    call(sf_rm, shell=True)
    p_tmp_rm="rm -rf picard_tmp"
    call(p_tmp_rm,shell=True)
    swift_cmd="swift upload PANCAN ./ --object-name ALIGN_TEST/" + bid + " -S " + str(ONE_GB)
    sys.stderr.write(date_time() + swift_cmd + "\n")
    call(swift_cmd,shell=True)
    cleanup_cmd="rm *.qs *bam* *.hist *.metrics *.txt *.bai"
    sys.stderr.write(date_time() + cleanup_cmd + "\n")
    call(cleanup_cmd, shell=True)

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
