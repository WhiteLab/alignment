#!/usr/bin/python
import sys
import json
import re
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
import subprocess

def parse_config(config_file):
    config_data=json.loads(open(config_file, 'r').read())
    return (config_data['tools']['samtools'],config_data['refs']['fa_ordered'])

def create_pos_ref(pos):
    fh=open(pos,'r')
    skip=next(fh)
    out=open('pos_list.txt','w')
    idx={}
    for line in fh:
        info=line.split('\t')
        if info[0] not in idx:
            idx[info[0]]={}
        idx[info[0]][info[1]]=1
    fh.close()
    for chrom in idx:
        for pos in idx[chrom]:
            out.write(chrom + '\t' + pos + '\n')
    out.close()

def base_call(cov,pile,base,samp,chrom,pos,i,k,m):
    if base not in cov[samp][chrom][pos]['bp']:
        cov[samp][chrom][pos]['bp'][base]={}
        cov[samp][chrom][pos]['bp'][base]['ct']=0
        cov[samp][chrom][pos]['bp'][base]['score']=[]
    cov[samp][chrom][pos]['bp'][base]['ct']+=1
    #convert ascii character of base score to phred and append
    cov[samp][chrom][pos]['bp'][base]['score'].append(ord(pile[i+1][m])-33)

def parse_pileup(out,tbl,sample,mode,cov,index):
    samp_list=[]
    # initilize index dictionary
    if mode=='b':
        fs=open(sample,'r')
        for line in fs:
            line=line.rstrip('\n')
            samp_list.append(line)
            index[line]={}
        fs.close()
    else:
        index[sample]={}
        samp_list.append(sample)
    fp=open(tbl,'r')
    next(fp)
    for line in fp:
        for samp in samp_list:
            line=line.rstrip('\n')
            (chrom,pos,ref,mut,gene)=line.split('\t')
            pos=int(pos)
            if chrom not in index[samp]:
                index[samp][chrom]={}
            if pos not in index[samp][chrom]:
                index[samp][chrom][pos]={}
            index[samp][chrom][pos]['ref']=ref
            # no list for ref, ought to be the same for all!
            index[samp][chrom][pos]['mut']=[]
            index[samp][chrom][pos]['mut'].append(mut)
            index[samp][chrom][pos]['gene']=gene

    fo=open(out,'r')    
    for line in fo:
        line=line.rstrip('\n')
        pile=line.split('\t')
        # fields 1-3 are chrom, pos, ref then triplets of ct, base, quality.
        (chrom,pos,ref)=(pile[0],int(pile[1]),pile[2])        
        k=0 # tracks sample list position
        for i in xrange(4,len(pile),3):
            if samp_list[k] not in cov:
                cov[samp_list[k]]={}
            if chrom not in cov[samp_list[k]]:
                cov[samp_list[k]][chrom]={}
            if pos not in cov[samp_list[k]][chrom]:
                cov[samp_list[k]][chrom][pos] ={}
                cov[samp_list[k]][chrom][pos]['bp']={}
            # marker to line up with quality scores to get per-base stats
            m=0 # tracks base quality score position within pile[i+1]
            # keep total coverage to check math
            cov[samp_list[k]][chrom][pos]['tot']=pile[i-1]
            j=0
            while j < len(pile[i]): # j is position of code for base call
                cur=pile[i][j] #current base
                # standard reference match
                call=0
                if re.match('[,|.]',cur):
                    base_call(cov,pile,ref,samp_list[k],chrom,pos,i,k,m)
                    call = 1
                    #convert ascii character of base score to phred and append
                # base substitution match
                elif re.match('[A|C|T|G|N|a|c|t|g|n]',cur):
                    cur=cur.upper()
                    base_call(cov,pile,cur,samp_list[k],chrom,pos,i,k,m)
                    call = 1
                elif cur == '^' or cur == '$':
                    # will make assumption that read will not start/end with in/del, skipping over read score
                    if cur =='^':
                        j=j+2
                    else:
                        j=j+1
                    cur=pile[i][j]
                    if re.match('[,|.]',cur):
                        base_call(cov,pile,ref,samp_list[k],chrom,pos,i,k,m)
                    else:
                        cur=cur.upper()
                        base_call(cov,pile,cur,samp_list[k],chrom,pos,i,k,m)
                    call = 1
                # making the final assumption that if none of the above conditions are met, then it's an indel
                elif cur == '-' or cur =='+':
                    size=pile[i][(j+1)]
                    try:
                        cur=pile[i][j] + pile[i][(j+2):(j+2+int(size))]
                    except:
                        sys.stderr.write('Parse error size was ' + str(size) + ' current base code ' + cur + ' chromosome ' + chrom + ' position ' + str(pos) + ' in sample ' + samp_list[k] + ' column index ' + str(i) + ' string index '+ str(j) + '\n')
                        exit(3)
                    j=j+2+int(size)
                    base_call(cov,pile,cur,samp_list[k],chrom,pos,i,k,m)
                    call = 1
                # odd * character seems to mean unaligned
                elif cur == '*':
                    cur='unaligned'
                    base_call(cov,pile,cur,samp_list[k],chrom,pos,i,k,m)
                    call = 1
                if call == 0:
                    sys.stderr.write('A base call was missed.  Check logic or input! current base code ' + cur + ' chromosome ' + chrom + ' position ' + str(pos) + ' in sample ' + samp_list[k] + ' column index ' + str(i) + ' string index '+ str(j) + '\n')
                    exit(3)
                m+=1
                j+=1
            k+=1
    return samp_list

# code obtained from web - may make it a class
def mean(data):
    """Return the sample arithmetic mean of data."""
    n = len(data)
    if n < 1:
        raise ValueError('mean requires at least one data point')
    return sum(data)/float(n) # in Python 2 use sum(data)/float(n)

def _ss(data):
    """Return sum of square deviations of sequence data."""
    c = mean(data)
    ss = sum((x-c)**2 for x in data)
    return ss

def pstdev(data):
    """Calculates the population standard deviation."""
    n = len(data)
    if n < 2:
        raise ValueError('variance requires at least two data points')
    ss = _ss(data)
    pvar = ss/n # the population variance
    return pvar**0.5

def calc_values(cov):
    for samp in sorted(cov):
        for chrom in sorted(cov[samp]):
            for pos in sorted(cov[samp][chrom]):
                tot_check=0
                for base in sorted(cov[samp][chrom][pos]['bp']):
                    score_tot=0
                    tot_check+=cov[samp][chrom][pos]['bp'][base]['ct']
                    #calc mean score
                    cov[samp][chrom][pos]['bp'][base]['x']="{0:.2f}".format(mean(cov[samp][chrom][pos]['bp'][base]['score']))
                    # calc std dev
                    if cov[samp][chrom][pos]['bp'][base]['ct'] >= 2:
                        cov[samp][chrom][pos]['bp'][base]['stdev']="{0:.2f}".format(pstdev(cov[samp][chrom][pos]['bp'][base]['score']))
                    else:
                        cov[samp][chrom][pos]['bp'][base]['stdev']='NA'
                if tot_check != int(cov[samp][chrom][pos]['tot']):
                    sys.stderr.write('Per base coverage total (' + str(tot_check) + ') does not equal total coverage ( ' + cov[samp][chrom][pos]['tot'] + ') of position (' + pos + ') for sample ' + samp + '.  Try using math.\n')
                    exit(69)
#                else:
#                    sys.stderr.write('Per base coverage total (' + str(tot_check) + ') equals total coverage ( ' + cov[samp][chrom][pos]['tot'] + ') of position (' + pos + ') for sample ' + samp + '.  Good job!\n')
def gen_report(cov,index,samp_list):
    sys.stdout.write('sample\tchromosome\tposition\tgene\treference\tbase\tct\tavg_score\tstd dev\ttracked_alt?\n')
    for samp in samp_list:
        for chrom in index[samp]:
            for pos in index[samp][chrom]:
                gene=index[samp][chrom][pos]['gene']
                ref=index[samp][chrom][pos]['ref']
                if pos in cov[samp][chrom]:
                    for base in cov[samp][chrom][pos]['bp']:
                        flag='N'
                        if base in index[samp][chrom][pos]['mut']:
                            flag = 'Y'
                        sys.stdout.write("\t".join((samp,chrom,str(pos),gene,ref,base,str(cov[samp][chrom][pos]['bp'][base]['ct']),str(cov[samp][chrom][pos]['bp'][base]['x']),str(cov[samp][chrom][pos]['bp'][base]['stdev']),flag)) + '\n')
                else:
                    sys.stdout.write("\t".join((samp,chrom,str(pos),gene,ref,'None',0,'NA','NA','NA')) + '\n')

def pre_report(mode,bam,sample,pos,config_file,ref_mnt):
    (samtools_tool,samtools_ref)=parse_config(config_file)
    samtools_ref = ref_mnt + '/' + samtools_ref
    create_pos_ref(pos)
    sys.stderr.write(date_time() + 'Creating mpileup with samtools\n')
    pre_rpt_cmd=samtools_tool + ' mpileup -D -d 500000 -l pos_list.txt -f ' +  samtools_ref
    out=''
    if mode == 'b':
        out='batch_pileup.txt'
        pre_rpt_cmd += ' -b  ' + bam +  ' > ' + out
    else:
        out=sample + '_pileup.txt'
        pre_rpt_cmd += ' ' + bam + ' > ' + out

    sys.stderr.write(date_time() + pre_rpt_cmd + "\n")
    try:
        subprocess.call(pre_rpt_cmd,shell=True)
    except:
        sys.stderr.write(date_time() + 'Pileup failed\n')
    cov={}
    index={}
    sys.stderr.write(date_time() + 'Parsing mpileup output\n')
    samp_list=parse_pileup(out,pos,sample,mode,cov,index)
    sys.stderr.write(date_time() + 'Calculating means and standard deviations of base quality scores\n')
    calc_values(cov)
    sys.stderr.write(date_time() + 'Generating report\n')
    gen_report(cov,index,samp_list)
    sys.stderr.write(date_time() + 'Report complete\n')
    return 0

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Variant pre-report using mpileup from samtools module.  Typically used with merged bams gives coverage report of mutation hot spots of interest.',add_help=True)
    parser.add_argument('-m','--mode',action='store',dest='mode',help='Run in (b)atch mode or (s)ingle sample mode')
    parser.add_argument('-b','--bam',action='store',dest='bam',help='List or name of merged bam file(s)')
    parser.add_argument('-s','--sample',action='store',dest='sample',help='Sample name or list')
    parser.add_argument('-p','--position',action='store',dest='pos',help='Tab-separated position list with header Chromosome position ref_base mut_base gene_sym')
    parser.add_argument('-j','--json',action='store',dest='config_file',help='JSON config file with tool and reference locations')
    parser.add_argument('-r','--ref_mount',action='store',dest='ref_mnt',help='Reference mount directory, i.e. /mnt/cinder/REFS_XXX')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (mode,bam,sample,pos,config_file,ref_mnt)=(inputs.mode,inputs.bam,inputs.sample,inputs.pos,inputs.config_file,inputs.ref_mnt)
    pre_report(mode,bam,sample,pos,config_file,ref_mnt)
