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