#!/usr/bin/env python2.7

## Obtained by MBROWN on git hun from https://github.com/cpwardell/bin/blob/master/metalfox.py, repurposed
## Purpose: Metal Fox is a variant filtering tool.

## Accepts either an exome directory as input, or
## full paths to the MuTect call_stats.out, vcf and bam files

## It calculates the following metrics with these limits:

# COMPLETE ## FoxoG must pass for C>A|G>T SNVs (see further details below) - using Pysam
# UNCERTAIN IF WORKING ## Minimum of 1 ALT read in each direction - using Pysam
# COMPLETE ## Mean Phred base quality must be > 26 - using Pysam
# COMPLETE ## Mean mapping quality must be >= 50 - using Pysam
# COMPLETE ## Alignability at site must be 1 - using Pysam

## Full FoxoG explanation:
## Use this with the MuTect tumor_lod and filter appropriately, as described in Costello et al 2013:
## "Discovery and characterization of artifactual mutations in deep coverage targeted capture
## sequencing data due to oxidative DNA damage during sample preparation"
## The cut-off is described is tumor_lod > -10+(100/3)FoxoG = OK
# DEVELOPING # We tweaked this for our data as tumor_lod > -19+(100/2)FoxoG = OK


import argparse
import logging
import os
import pysam
import csv
import sys
import pybedtools

## Gather command line args
parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, help="full path to exome directory",required=False)
parser.add_argument("-f1", type=str, help="full path to MuTect call_stats.out file",required=False)
parser.add_argument("-f3", type=str, help="full path to bam",required=False)
parser.add_argument("-m", type=str, help="mapability file",required=False)

parser.add_argument("--debug", help="Turn debugging mode on",action="store_true")
args = parser.parse_args()

## Turn logging on if desired
if(args.debug):
    logging.basicConfig(level=logging.DEBUG)
    logging.debug("Debugging mode enabled")

## If args.f is supplied, infer args.f1,args.f3
if(args.f!=None):
    args.f1 = os.path.join(args.f,"mutect/call_stats.out")
    args.f3 = os.path.join(args.f,"dedup/dedup.bam")

## If args.f1 and args.f3 don't exist by this point, they haven't been explicitly supplied
## or inferred, so warn the user and exit
if(args.f1==None and args.f3==None):
    print "You must supply EITHER the -f OR -f1, -f3 arguments"
    sys.exit()

## Check files exist.  If they don't, warn user and exit
if not os.path.isfile(args.f1):
    print args.f1+" does not exist"
    sys.exit()
if not os.path.isfile(args.f3):
    print args.f3+" does not exist"
    sys.exit()

## Open the bamfile.  We do it ONCE and ONCE ONLY for speed reasons
logging.debug("Opening bam file")
samfile=pysam.Samfile(args.f3,"rb") # rb = "read bam"

## Define a list in which to store the indices of all the lines we want to keep
## Currently unused, but may be useful later...
logging.debug("Calculating foxog")
foxkeepers = []
logging.debug("Calculating ALT read directions")
directionkeepers=[]
logging.debug("Calculating map and base quality scores")
mapbasekeepers=[]
logging.debug("Calculating alignability")
alignkeepers=[]

## Function for calculating the mean; we don't want to use numpy
def mean(numbers):
    x=float(sum(numbers))/len(numbers)
    return(x)

## Iterate through every SNV
## We use enumerate() to create a nice index for us to use
for idx,row in enumerate(csv.reader(open(args.f1),delimiter="\t")):
    ## Output row number for debugging
    if(idx%1000==0):
	logging.debug(str(idx)+" rows analysed")

    ## Skip the first row, print the header row and skip rows with REJECT in them
    if(row[0].startswith("#")):
	continue
    if(row[0].startswith("contig")):
	header="\t".join(row)
	print "\t".join([header,"FoxoG","forward_reads","reverse_reads","mean_mapq",
	"mean_baseq","alignability"])
	continue
    if("KEEP" not in row[50]):
	continue

    ## Set properties of SNV
    CHROM=row[0]
    POS=int(row[1])-1 # Pysam coordinates are ZERO-based, so we MUST subtract 1
    REF=row[3]
    ALT=row[4]
    T_LOD_FSTAR=float(row[20])

    ## Get values at the SNV location
    F2R1=float() # Orientation for FoxoG calc
    F1R2=float() # Orientation for FoxoG calc
    F1=0 # Orientation for orientation filter
    F2=0 # Orientation for orientation filter
    mapq=[] # mapping quality scores
    baseq=[] # base quality scores
    FoxoG="NA"
    alignability="NA"

    for alignedread in samfile.fetch(CHROM,POS,POS+1):
        if not alignedread.is_proper_pair:
	    continue
	else:
	    ## Which base in the read is at the position we want?  Use the
	    ## "aligned_pairs" list of tuples to determine this
	    offset = [item for item in alignedread.aligned_pairs if item[1] == POS][0][0]
	    if(offset!=None):
		mapq.append(alignedread.mapq)
		## We subtract 33 because SAM specification tells us to
		baseq.append(ord(alignedread.qual[offset])-33)
	    if(offset!=None and alignedread.seq[offset]==ALT):
		if(alignedread.is_read1 and alignedread.is_reverse): # 83/pPr1 F2R1
		    F2R1+=1
		if(alignedread.is_read2 and alignedread.mate_is_reverse): # 163/pPR2 F2R1
		    F2R1+=1
		if(alignedread.is_read1 and alignedread.mate_is_reverse): # 99/pPR1 F1R2
		    F1R2+=1
		if(alignedread.is_read2 and alignedread.is_reverse): # 147/pPr1 F1R2
		    F1R2+=1
		if(alignedread.is_reverse):
		    F1+=1
		if(alignedread.mate_is_reverse):
		    F2+=1
    ## We skip positions with no correctly paired ALT reads
    if(len(mapq)==0 or len(baseq)==0):
	continue
    elif(mean(mapq)>=50 and mean(baseq)>=25):
	mapbasekeepers.append(idx)
    ## Calculation of FoxoG - only necessary if C>A|G>T SNV
    ## Equation is: ALT_F1R2/(ALT_F1R2 + ALT_F2R1) or ALT_F2R1/(ALT_F1R2 + ALT_F2R1)
    ## C>anything:  numerator is ALT_F2R1
    ## A>anything:  numerator is ALT_F2R1
    ## G>anything:  numerator is ALT_F1R2
    ## T>anything:  numerator is ALT_F1R2
    if((row[3]=="C" and row[4]=="A") or (row[3]=="G" and row[4]=="T")):
	## If sum of F1R2 and F2R1 is zero, all reads have an indel in them, so it should be removed
	if((F1R2 + F2R1)!=0):
	    if(REF=="C"):
		FoxoG = F2R1/(F1R2 + F2R1)
	    if(REF=="G"):
		FoxoG = F1R2/(F1R2 + F2R1)
	## If FoxoG is still "NA" at this point, it must be rubbish, so set it to 1
	if(FoxoG=="NA"):
	    FoxoG=1
	#score=-19+(100/2)*FoxoG
	score=-10+(100/3)*FoxoG
	#logging.debug("\t".join([str(FoxoG),str(T_LOD_FSTAR),str(score)]))
	if(T_LOD_FSTAR>score):
	    foxkeepers.append(idx)
    else:
	foxkeepers.append(idx)

    ## Calculate directionkeepers metric
    if(F1>=1 and F2>=1):
	directionkeepers.append(idx)

    ## We put this in a try clause, as some chromosome names (MT, GL, etc) are
    ## not in the alignability file and generate errors

    alignfile=pysam.Tabixfile(args.m)
    try:
	for record in alignfile.fetch(CHROM, POS, POS+1):
	    alignability=float(record.split("\t")[3])
	    if(alignability>=1): # is alignability 1?
		alignkeepers.append(idx)
    except:
	pass

    ## OUTPUT RESULTS HERE
    tabrow = "\t".join(row)
    print "\t".join([tabrow,str(FoxoG),str(F1),str(F2),str(mean(mapq)),str(mean(baseq)),str(alignability)])

## Close the bamfile
samfile.close()