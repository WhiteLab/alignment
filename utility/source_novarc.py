#!/usr/bin/env python
import os
import sys


def source_novarc(novarc):
    with open(novarc, 'r') as f:
        for line in f:
            k, v = line.rstrip().split('=')
            k = k.replace('export ', '')
            v = v.replace('"', '')
            os.environ[k] = v


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Standalone method to source openstack variables when interacting with swift')
    parser.add_argument('-n', '--novarc', action='store', dest='novarc', help='.novarc with openstack authentication')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (novarc) = (inputs.novarc)
    source_novarc(novarc)