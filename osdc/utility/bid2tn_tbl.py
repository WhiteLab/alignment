#!/usr/bin/python

import sys
fh=open(sys.argv[1],'r')

for line in fh:
    tumor=line.rstrip('\n')
    titems=line.split('\t')
    normal=next(fh)
    normal=normal.rstrip('\n')
    nitems=normal.split('\t')
    sys.stdout.write('\t'.join((titems[0] + '_' + nitems[0],titems[0],nitems[0])) + '\n')
fh.close()
