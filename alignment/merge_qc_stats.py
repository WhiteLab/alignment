#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')


def download_from_swift(p_dir, f_dir, lane_list):
    lanes = open(lane_list, 'r')
    head = ''
    data = []
    for line in lanes:
        line = line.rstrip('\n')
        (bid, seqtype, lane_csv) = line.split('\t')
        for lane in lane_csv.split(', '):
            cur = p_dir + '/' + f_dir + '/' + bid + '/QC/' + bid + '_' + lane + '.qc_stats.txt'
            try:
                stat = open(cur, 'r')
                head = next(stat)
                data.append(next(stat))
                stat.close()
            except:
                sys.stderr.write('Could not find file ' + cur + ' skipping!\n')
    lanes.close()
    sys.stdout.write(head)
    for datum in data:
        sys.stdout.write(datum)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Uses pipeline lane list to create a summary table of qc stats')
    parser.add_argument('-p', '--project_dir', action='store', dest='p_dir', help='Project directory, i.e. '
                                                                                  '/cephfs/PROJECTS/PANCAN')
    parser.add_argument('-d', '--file_dir', action='store', dest='f_dir',
                        help='file directory prefix with qc files')
    parser.add_argument('-l', '--lane_list', action='store', dest='lane_list',
                        help='Original lane list used to run pipeline')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (p_dir, f_dir, lane_list) = (inputs.p_dir, inputs.f_dir, inputs.lane_list)
    download_from_swift(p_dir, f_dir, lane_list)
