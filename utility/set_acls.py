#!/usr/bin/env python3
import sys
import subprocess


def set_acls(dirname, user, group):
    cmd = 'chown -R ' + user + ':' + group + ' ' + dirname
    sys.stderr.write('Setting acls ' + cmd + '\n')
    subprocess.call(cmd, shell=True)
