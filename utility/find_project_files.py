#!/usr/bin/env python

import subprocess
import sys
from subprocess import check_output

from date_time import date_time


def find_project_files(project, subdir):
    find_cmd = "find " + project + "/" + subdir
    sys.stderr.write(date_time() + find_cmd + "\n")
    try:
        results = check_output(find_cmd, shell=True, stderr=subprocess.PIPE)
        return results
    except:
        sys.stderr.write(date_time() + "Search of " + subdir + " from " + project + " failed\n")
        exit(1)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Simple download module to get files from swift.  Can use prefix or whole object name')
    parser.add_argument('-p', '--project', action='store', dest='project', help='Project directory')
    parser.add_argument('-s', '--sub-directory', action='store', dest='subdir',
                        help='File subdirectory within project')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (project, subdir) = (inputs.project, inputs.subdir)
    find_project_files(project, subdir)
