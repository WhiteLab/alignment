#!/usr/bin/python
import sys
#import pdb
hist={}

i=1
m=1000000
score_type=sys.argv[1]
adj={}
# Sanger, Solexa, Illumina 1.3+, Illumina 1.5+, Illumina 1.8+
adj['S']=-33
adj['X']=-64
adj['I']=-64
adj['J']=-64
adj['L']=-64
for line in sys.stdin:
    if i % m == 0:
        sys.stderr.write('At line ' + str(i) + '\n')
    line=line.rstrip('\n')
    fields=line.split()

    qual=fields[10]
    for base in qual:
        score=(ord(base) + adj[score_type])
        if score not in hist:
            hist[score]=1
        else:
            hist[score]+=1
    i+=1
for score in sorted(hist):
    sys.stdout.write(str(score) + '\t' + str(hist[score]) + '\n')
