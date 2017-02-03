#!/usr/bin/env python

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from analysis.variant_annot_pipe import *


def override_rejected_variants(config_file, table, ref_mnt):
    var_dict = {}
    # hash variants to recover
    for line in open(table):
        # format pair\tchrom\tpos\tref\talt\ttype
        (pair, chrom, pos, ref, alt, ty) = line.rstrip('\n').split('\t')
        if pair not in var_dict:
            var_dict[pair] = {}
        ty = 'snv'
        if len(ref) > 1 or len(alt) > 1:
            ty = 'indel'
        if ty not in var_dict[pair]:
            var_dict[pair][ty] = []
        var_dict[pair][ty].append('\t'.join((chrom, pos, ref, alt)))



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Recover rejected variants deemed to be real by manual curation and '
                                                 'reincoprorate into reports.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-t', '--table', action='store', dest='table', help='Table with pairs and coords to recover')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference mount directory, i.e. /mnt/cinder/REFS_XXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, table, ref_mnt) = (inputs.config_file, inputs.table, inputs.ref_mnt)
    override_rejected_variants(config_file, table, ref_mnt)
