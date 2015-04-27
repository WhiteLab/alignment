#!/usr/bin/python
import sys
import re
from date_time import date_time
from subprocess import call
from subprocess import check_output
import subprocess 
import glob

def upload_variants_to_swift(cont,obj,sample_list,sample_pairs):
   src_cmd='. ~/.novarc;'
   ONE_GB = 1073741824
   fh=open(sample_list,'r')
   for sample in fh:
      sample = sample.rstrip('\n')
      swift_cmd=src_cmd + 'swift upload ' + cont + ' BAM/' + sample + '.merged.final.bam -S ' + str(ONE_GB) + ' --skip-identical --object-name ' + obj + '/' + sample + '/BAM/' + sample + '.merged.final.bam >> LOGS/' + sample + '.upload.log 2>> LOGS/' + sample + '.upload.log'
      check=call(swift_cmd,shell=True)
      swift_cmd=src_cmd + 'swift upload ' + cont + ' BAM/' + sample + '.merged.final.bai -S ' + str(ONE_GB) + ' --skip-identical --object-name ' + obj + '/' + sample + '/BAM/' + sample + '.merged.final.bai >> LOGS/' + sample+ '.upload.log 2>> LOGS/' + sample + '.upload.log'
      try:
         call(swift_cmd,shell=True)
      except:
         swift_cmd=src_cmd + 'swift upload ' + cont + ' BAM/' + sample + '.merged.final.bam.bai -S ' + str(ONE_GB) + ' --skip-identical --object-name ' + obj + '/' + sample + '/BAM/' + sample + '.merged.final.bai >> LOGS/' + sample+ '.upload.log 2>> LOGS/' + sample + '.upload.log'
         check=check + call(swift_cmd,shell=True)
      if check==0:
         sys.stderr.write(date_time() + 'Uploading final BAMs for ' + sample + ' successful!\n')
      else:
         sys.stderr.write(date_time() + 'Uploading final BAMs for ' + sample + ' failed\n')
         exit(1)
   fh.close()

   # upload analysis and annotation files                                                                                                                                                                                         
   suffix_list1=['.out','.out.keep','.vcf','.vcf.keep']
   suffix_list2=['.out.keep.eff.vcf','.vcf.keep.eff.xls']
   fh=open(sample_pairs,'r')
   for line in fh:
      pair = line.split('\t')[0]
      # upload analysis files
      for suffix in suffix_list1:
         swift_cmd=src_cmd + 'swift upload ' + cont + ' ANALYSIS/' + pair + suffix + ' -S ' + str(ONE_GB) + ' --skip-identical --object-name ANALYSIS/' + pair + '/OUTPUT/' + pair + suffix + ' >> LOGS/' + pair + '.upload.log 2>> LOGS/' + pair + '.upload.log'
         check=call(swift_cmd,shell=True)
         if check==0:
            sys.stderr.write(date_time() + 'Uploading analysis file ' + pair + suffix + ' successful!\n')
         else:
            sys.stderr.write(date_time() + 'Uploading analysis file ' + pair + suffix + ' failed!\n')
            exit(1)
      # upload annotation files
      for suffix in suffix_list2:
         swift_cmd=src_cmd + 'swift upload ' + cont + ' ANNOTATION/' + pair + suffix + ' -S ' + str(ONE_GB) + ' --skip-identical --object-name ANNOTATION/' + pair + '/OUTPUT/' + pair + suffix + ' >> LOGS/' + pair + '.upload.log 2>> LOGS/' + pair + '.upload.log'
         check=call(swift_cmd,shell=True)
         if check==0:
            sys.stderr.write(date_time() + 'Uploading annotation file ' + pair + suffix + ' successful!\n')
         else:
            sys.stderr.write(date_time() + 'Uploading annotation file ' + pair + suffix + ' failed!\n')
            exit(1)
      # upload log files
      mut_list=glob.glob('LOGS/' + pair + '.mut*')
      for mut in mut_list:
         swift_cmd=src_cmd + 'swift upload ' + cont + ' ' + mut + ' -S ' + str(ONE_GB) + ' --skip-identical --object-name ANALYSIS/' + pair + '/' + mut +' >> LOGS/' + pair + '.upload.log 2>> LOGS/' + pair + '.upload.log'
         check=call(swift_cmd,shell=True)
         if check==0:
            sys.stderr.write(date_time() + 'Uploading analysis log file ' + mut  + ' successful!\n')
         else:
            sys.stderr.write(date_time() + 'Uploading analysis log file ' + mut + ' failed!\n')
            exit(1)
      swift_cmd=src_cmd + 'swift upload ' + cont + ' LOGS/' + pair + '.snpeff.log -S ' + str(ONE_GB) + ' --skip-identical --object-name ANNOTATION/' + pair + '/LOGS/' + pair + '.snpeff.log >> LOGS/' + pair + '.upload.log 2>> LOGS/' + pair + '.upload.log'
      check=call(swift_cmd,shell=True)
      if check==0:
         sys.stderr.write(date_time() + 'Uploading annotation log file ' + pair + '.snpeff.log' + ' successful!\n')
      else:
         sys.stderr.write(date_time() + 'Uploading annotation log file ' + pair + '.snpeff.log' + ' failed!\n')
         exit(1)
   fh.close()
    

   return 0
if __name__ == "__main__":
   import argparse
   parser=argparse.ArgumentParser(description='Uploads current directory contents to specified object and container')
   parser.add_argument('-c','--container',action='store',dest='cont',help='Swfit container name to upload to.  i.e. PANCAN')
   parser.add_argument('-o','--object',action='store',dest='obj',help='Swift object name root to use for aligned merged bam files.  i.e. ALIGN/2015-1234')
   parser.add_argument('-sl','--sample_list',action='store',dest='sample_list',help='Sample list, one per line')
   parser.add_argument('-sp','--sample_pairs',action='store',dest='sample_pairs',help='Sample tumor/normal pairs, tsv file with bid pair, sample1, sample2')
   
   if len(sys.argv)==1:
      parser.print_help()
      sys.exit(1)
      
   inputs=parser.parse_args()
   (cont,obj,sample_list,sample_pairs)=(inputs.cont,inputs.obj,inputs.sample_list,inputs.sample_pairs)
   upload_variants_to_swift(cont,obj,sample_list,sample_pairs)
