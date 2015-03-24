#!/usr/bin/python

import sys
import pdb
import re
import subprocess

#curl -X PUT -H 'Content-Type: application:json' https://128.135.219.167:6984/igsb_qc/24ec4b43cfe304ff4709e76f74004da8?rev=1-48f57c4e6aaba987d0786e6750c19d7c -d @ALIGN_TEST/2014-2241/QC/2014-2241_141126_SN1070_0304_BHB5BPADXX_2.qc_stats.json -k -u "mbrown:Jah6Eeji"
#ALIGN_TEST/2014-2240/QC/2014-2240_141126_SN1070_0304_BHB5BPADXX_2.qc_stats.json{"ok":true,"id":"24ec4b43cfe304ff4709e76f7400074d","rev":"1-7dd84750ce781a128cfdb7e14ae4d773"}
curl='curl -X PUT -H "Content-Type: application:json" https://128.135.219.167:6984/igsb_qc/'
login=' -k -u "mbrown:Jah6Eeji"'
fh=open(sys.argv[1],'r')
setproxy='export http_proxy=http://cloud-proxy:3128; export https_proxy=http://cloud-proxy:3128; export no_proxy="rados-bionimbus-pdc.opensciencedatacloud.org";'
for line in fh:
#    pdb.set_trace()
    fields=line.split('\t')
    f=re.search('"id":"(.*)","rev":"(.*)"',fields[1])
    c_cmd=setproxy + curl + f.group(1) + '?rev=' + f.group(2) +  ' -d @' + fields[0] + login
    sys.stderr.write(c_cmd + '\n')
    check=subprocess.check_output(c_cmd,shell=True)
    sys.stdout.write(fields[0] + '\t' + check)
