def create_index(ref):
    # helper script common to report generators
    fh = open(ref)
    g_dict = {}
    for line in fh:
        info = line.rstrip('\n').split('\t')
        g_dict[info[-1]] = info[-2]
        g_dict[info[1]] = info[-2]
    fh.close()
    return g_dict


def create_target(cfile):
    import re
    cindex = {}
    fh = open(cfile, 'r')
    for line in fh:
        m = re.search('(chr\w+):(\d+)-(\d+)', line)
        try:
            (chrom, start, end) = (m.group(1), m.group(2), m.group(3))
        except:
            sys.stderr.write(line + ' doesn\'t fit format (chr\w+):(\d+)-(\d+), skipping\n')
            continue
        if chrom not in cindex:
            cindex[chrom] = {}
        cindex[chrom][(int(start))] = int(end)
    return cindex


def mark_target(chrom, pos, on_dict):
    f = 0
    if chrom in on_dict:
        for start in sorted(on_dict[chrom]):
            if start <= int(pos) <= on_dict[chrom][start]:
                f = 1
                break
            elif start > int(pos):
                break
    status = "OFF"
    if f == 1:
        status = "ON"
    return status


def calc_pct(a, b):
    # return both formatted and unformatted
    ratio = float(b) / (float(a) + float(b)) * 100
    fmt = "{0:.2f}%".format(ratio)
    return ratio, fmt

