#!/usr/bin/perl -w
use strict;
use warnings;

my $pair_list = shift or die "Usage: <file list of pairs, 3 columns, sample - idT - idN\n";

open PAIRS, "<$pair_list";
my $home = "/raid/users/jgrundst";
my $ref = "$home/REF";
my $ref_fa = "$ref/hg19_Ordered.fa";
#my $dbsnp = "$ref/GATK_hg19_bundle/dbsnp_137.hg19.ordered.vcf";
#my $cosmic = "$ref/COSMIC/Cosmic.hg19.v2.vcf";
#my $intervals_to_process = "$ref/nimblegen.merged.intervals";
my $intervals_to_process = "$ref/hg19_Ordered.intervals";
#my $java = "$home/jdk1.7.0_45/bin/java";
my $java = "/usr/bin/java";
my $tmp = "-Djava.io.tmpdir=$home/temp";
my $run_muTect = "$java $tmp -Xmx2g -jar $home/TOOLS/muTect-1.1.4.jar";
`mkdir -p scripts/`;
`mkdir -p logs/`;
#list all the sample names, and where the sample number matches i.e. 201T and 201B, output the bamfiles, one for blood and on for tumor so they can go into mpileup and bcftools
foreach my $line (<PAIRS>){
	chomp $line;
	(my $s, my $tumor_id, my $normal_id) = split /\t/, $line;
	my $normal_bam = `ls $normal_id*.rmdup.srt.bam`;
	my $tumor_bam = `ls $tumor_id*.rmdup.srt.bam`;
	chomp $normal_bam;
	chomp $tumor_bam;

	print STDERR "T: $tumor_bam	N: $normal_bam\n";
	my $out = "$tumor_id-$normal_id";
	foreach my $interval (`cat $intervals_to_process`) {
		chomp $interval;
		print "interval: ->$interval<-\n";
		(my $chr, my $garbage) = split /:/, $interval;
		my $output_file = "$out.$chr.out";
		my $vcf_file = "$out.$chr.vcf";
		my $log_file = "$out.$chr.log";
		my $run_file = "scripts/run_$out.$chr.sh";

		`echo "#!/bin/bash
#PBS -N $out.$chr
#PBS -l nodes=1:ppn=2
#PBS -j oe
#PBS -o logs/
cd \\\$PBS_O_WORKDIR
$run_muTect -T MuTect -R $ref_fa --intervals $interval --input_file:normal $normal_bam --input_file:tumor $tumor_bam --out $output_file -vcf $vcf_file --enable_extended_output --strand_artifact_power_threshold 0 -log $log_file
cat $output_file | grep -v REJECT > $output_file.keep
cat $vcf_file | grep -v REJECT > $vcf_file.keep" > $run_file`;
		`chmod +x $run_file`;
	}

}
