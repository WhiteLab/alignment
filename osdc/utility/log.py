import sys

def log(loc,string):
    lh=open(loc,'a')
    lh.write(string)
    lh.flush()
if __name__ == "__main__":
    import argparse

    parser=argparse.ArgumentParser(description='Tiny tool to direct log output of modules.  Really meant to be used within and not on its own')
    parser.add_argument('-l','--location',action='store',dest='loc',help='Location to write output to')
    parser.add_argument('-s','--string',action='store',dest='string',help='String to output')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    loc=inputs.loc
    string=inputs.string
    log(loc,string)
