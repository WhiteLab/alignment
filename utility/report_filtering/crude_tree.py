#!/usr/bin/env python
import sys
import pdb


def crude_tree(stable, vtable):
    sdict = {}
    gdict = {}
    samples = open(stable)
    head = next(samples)
    group = head.rstrip('\n').split('\t')
    #pdb.set_trace()
    for pairs in samples:
        pset = pairs.rstrip('\n').split('\t')
        for i in xrange(0, len(pset), 1):
            if len(pset[i]) > 1:
                sdict[pset[i]] = group[i]
                if group[i] not in gdict:
                    gdict[group[i]] = {}
                    gdict[group[i]]['ct'] = 0
                gdict[group[i]]['ct'] += 1
    samples.close()
    vars = open(vtable)
    samps = next(vars)
    slist = samps.rstrip('\n').split('\t')
    vlist = []
    for var in vars:
        vals = var.rstrip('\n').split('\t')
        cur_var = vals[0]
        vlist.append(vals[0])
        for i in xrange(1, len(vals), 1):
            if float(vals[i]) > 0:
                cur_gr = sdict[slist[i]]
                if cur_var not in gdict[cur_gr]:
                    gdict[cur_gr][cur_var] = 0
                gdict[cur_gr][cur_var] += 1
    vars.close()
    sys.stdout.write('Group/Variant Count')
    for gr in group:
        sys.stdout.write('\t' + gr + '_n=' + str(gdict[gr]['ct']))
    print
    for var in vlist:
        sys.stdout.write(var)
        for gr in group:
            if var in gdict[gr]:
                sys.stdout.write('\t' + str(gdict[gr][var]))
            else:
                sys.stdout.write('\t0')
        print


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Creates table with variant occurrence counts by sample group.')
    parser.add_argument('-t', '--table', action='store', dest='stable',
                        help='Sample table')
    parser.add_argument('-v', '--variant', action='store', dest='vtable',
                        help='Variant table')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (stable, vtable) = (inputs.stable, inputs.vtable)
    crude_tree(stable, vtable)