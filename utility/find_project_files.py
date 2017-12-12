#!/usr/bin/env python3

import subprocess
import sys
from subprocess import check_output
import pdb
from utility.date_time import date_time


def find_project_files(file_dir, file_prefix):
    find_cmd = "find " + file_dir + " -name " + file_prefix + '*'
    sys.stderr.write(date_time() + find_cmd + "\n")
    try:
        results = check_output(find_cmd, shell=True, stderr=subprocess.PIPE)
        pdb.set_trace()
        return results
    except:
        sys.stderr.write(date_time() + "Search of " + file_prefix + " from " + file_dir + " failed\n")
        exit(1)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Simple download module to get files from swift.  Can use prefix or whole object name')
    parser.add_argument('-d', '--file-dir', action='store', dest='file_dir', help='Directory with relevant files')
    parser.add_argument('-p', '--file-prefix', action='store', dest='file_prefix',
                        help='Prefix of relevant files')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (file_dir, file_prefix) = (inputs.file_dir, inputs.file_prefix)
    find_project_files(file_dir, file_prefix)
