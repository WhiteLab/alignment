#!/usr/bin/python
# written by Miguel Brown 2015-Feb-11. test run based on align.pl calls just to confirm successful installation of tools on pipeline vm

import sys
import re
from subprocess import call

end1=sys.argv[1]
end2=sys.argv[2]

#inputs

SAMPLES={}

s=re.match('^(\S+)_1_sequence\.txt\.gz$',end1)

sample=s.group(1)
HGACID=sample.split("_")
SAMPLES[sample]={}
SAMPLES[sample]['f1']=end1
SAMPLES[sample]['f2']=end2
RGRP="@RG\\tID:" + sample + "\\tLB:" + HGACID[0] + "\\tSM:" + HGACID[0] + "\\tPL:illumina"
#RGRP="@RG\tID:" + sample + "\tLB:" + HGACID[0] + "\tSM:" + HGACID[0] + "\tPL:illumina"
print RGRP;

#tools and refs

fastx='fastx_quality_stats'
bwa='/home/ubuntu/TOOLS/bwa-0.7.8/bwa'
bwa_ref='/mnt/cinder/REFS/bwa-0.7.8/hg19.fa'
samtools='/home/ubuntu/TOOLS/samtools-0.1.19/samtools'
samtools_ref='/mnt/cinder/REFS/samtools-0.1.19/hg19.fa'
java='/home/ubuntu/TOOLS/jdk1.7.0_45/bin/java'
picard='/home/ubuntu/TOOLS/picard/dist/picard.jar'
picard_tmp='picard_tmp'
bedtools='/home/ubuntu/TOOLS/bedtools2/bin/bedtools'
exome_bed_ref='/mnt/cinder/REFS/BED/refseq.Hs19.coding_regions.merged.bed'
genome_bed_ref='/mnt/cinder/REFS/BED/hg19_complete_sorted.bed'
capture_bed_ref='/mnt/cinder/REFS/BED/capture_panel_2.0.bed'

fastx_cmd='gzip -dc ' + end1 + ' | ' + fastx + ' -N -o ' + sample + '_1.qs 2> fq1_log.txt'
#call(fastx_cmd,shell=True);
fastx_cmd='gzip -dc ' + end2 + ' | ' + fastx + ' -N -o ' + sample + '_2.qs 2> fq2_log.txt'
#call(fastx_cmd,shell=True);

bwa_cmd="(" + bwa + " mem -t 8 -R \"" + RGRP + "\" -v 2 " + bwa_ref + " " + end1 + " " + end2 + " | " + samtools + " view -bT " + samtools_ref + " - > " + sample + ".bam) > " + sample + ".bwa.pe.log 2>&1"

#call(bwa_cmd,shell=True)
print bwa_cmd

picard_sort_pe_cmd=java + " -Xmx8g -jar " + picard + " SortSam CREATE_INDEX=true TMP_DIR=" + picard_tmp + " INPUT=" + sample + ".bam OUTPUT=" + sample + ".srt.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT > " + sample + ".picard.sort.pe.log 2>&1"

print picard_sort_pe_cmd

picard_rmdup_cmd=java + " -Xmx8g -jar " + picard + " MarkDuplicates CREATE_INDEX=true TMP_DIR=" + picard_tmp + " REMOVE_DUPLICATES=true ASSUME_SORTED=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500 INPUT=" + sample + ".srt.bam OUTPUT=" + sample + ".rmdup.srt.bam METRICS_FILE=" + sample + ".rmdup.srt.metrics VALIDATION_STRINGENCY=LENIENT > " + sample + ".picard.rmdup.pe.log 2>&1"
print picard_rmdup_cmd

flagstats_cmd=samtools + " flagstat " + sample + ".srt.bam > " + sample + ".srt.bam.flagstats"
print flagstats_cmd

flagstats_cmd=samtools + " flagstat " + sample + ".rmdup.srt.bam > " + sample + ".rmdup.srt.bam.flagstats"
print flagstats_cmd

exome_coverage_cmd=bedtools + " coverage -abam " + sample + ".rmdup.srt.bam -b " + exome_bed_ref + " | grep all > " + sample + ".exome.hist"
print exome_coverage_cmd

genome_coverage_cmd=bedtools + " genomecov -ibam " + sample + ".rmdup.srt.bam -g " + genome_bed_ref + " | grep genome > " + sample + ".genome.hist"
print genome_coverage_cmd

target_coverage_cmd=bedtools + " coverage -hist -abam " + sample + ".rmdup.srt.bam -b " + capture_bed_ref + " | grep all > " + sample + ".capture.hist"
print target_coverage_cmd
