#!/usr/bin/env python3
import json
import os
import re
import sys
from subprocess import check_output

from date_time import date_time


def set_proxy(http_proxy, https_proxy, no_proxy):
    os.environ.update({'http_proxy': http_proxy, 'https_proxy': https_proxy, 'no_proxy': no_proxy})


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['server'], config_data['user'], config_data['password'], config_data['db'], \
           config_data['http_proxy'], config_data['https_proxy'], config_data['no_proxy']


def update_couchdb(fn, config_file):
    (server, user, password, db, http_proxy, https_proxy, no_proxy) = parse_config(config_file)
    set_proxy(http_proxy, https_proxy, no_proxy)
    fh = open(fn, 'r')
    # Get uuid command

    for obj in fh:
        obj = obj.rstrip('\n')
        get_uuid = 'curl -X GET ' + server + '/_uuids -k;'
        sys.stderr.write(date_time() + get_uuid + '\n')
        uuid_out = check_output(get_uuid, shell=True)
        m = re.findall('\"(\w+)\"', uuid_out)
        uuid = m[1]
        # typical response: {"uuids":["24ec4b43cfe304ff4709e76f7400074d"]}
        curl = 'curl -X PUT -d @' + obj + ' "' + server + '/' + db + '/' + uuid \
               + '" -H "Content-Type: application/json" -k -u "' + user + ':' + password + '"'
        couch_cmd = curl
        # get response
        sys.stderr.write(date_time() + couch_cmd + '\n')
        result = check_output(couch_cmd, shell=True)
        result = result.rstrip('\n')
        sys.stderr.write(obj + '\t' + result + '\n')
        if result == 1:
            sys.stderr.write(date_time() + 'Database update failed for qc stats.  Check connection')
            exit(1)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Update couch db with qc stats using a json object list file')
    parser.add_argument('-f', '--file', action='store', dest='fn', help='qc_stats.json document list')
    parser.add_argument('-j', '--json', action='store', dest='config_file', help='Config file with server info')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (fn, config_file) = (inputs.fn, inputs.config_file)
    update_couchdb(fn, config_file)
