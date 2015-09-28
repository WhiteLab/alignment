#!/usr/bin/perl -w
use strict;

my $SAMTOOLS='/raid/users/jgrundst/TOOLS/samtools-1.1/samtools';
my $VARSCAN='/glusterfs/users/JGRUNDSTAD/TOOLS/VarScan.v2.3.6.jar';
my $REF='/glusterfs/users/JGRUNDSTAD/REF/hg19/hg19_Ordered.fa';
my $INTERVALS='/glusterfs/users/JGRUNDSTAD/REF/hg19/hg19_Ordered.intervals';

# get bam list
my $SAMPLE_LIST=join ' ', @ARGV;

# create sample.txt for VarScan
if(-e 'samples.txt') { `rm samples.txt`; }
foreach my $s (@ARGV) {
	if($s =~ /^(2\d+-\d+)(\.|-|_)/) {
		`echo "$1" >> samples.txt`;
	}
}

`mkdir -p logs`;
my $out_dir = 'germline_vcfs';
`mkdir -p $out_dir`;
my $scripts = 'germline_scripts';
`mkdir -p $scripts`;

my $HEADER = "#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -j oe
#PBS -o logs/";

# Run VarScan mpileup2snp to call SNVs:
foreach my $region (`cat $INTERVALS`) {
	print "$region";
	chomp $region;
	(my $chr, my $range) = split /:/, $region;
	open SNP, ">$scripts/snp.$chr.sh" or die "$!";
	print SNP "$HEADER
#PBS -N snp.$chr
cd \$PBS_O_WORKDIR
$SAMTOOLS mpileup -B -q 1 -f $REF -r \"$region\" $SAMPLE_LIST | java -Xmx2g -jar $VARSCAN mpileup2snp --vcf-sample-list samples.txt --min-coverage 10 --min-var-freq 0.10 --p-value 0.05 --output-vcf 1 > $out_dir/snp.$chr.vcf 2> logs/mpileup2snp.$chr.log
";
	open INDEL, ">$scripts/indel.$chr.sh" or die "$!";
	print INDEL "$HEADER
#PBS -N indel.$chr
cd \$PBS_O_WORKDIR
$SAMTOOLS mpileup -B -q 1 -f $REF -r \"$region\" $SAMPLE_LIST | java -Xmx2g -jar $VARSCAN mpileup2indel --vcf-sample-list samples.txt --min-coverage 10 --min-var-freq 0.10 --p-value 0.10 --output-vcf 1 > $out_dir/indel.$chr.vcf 2> logs/mpileup2indel.$chr.log
";
	close SNP; close INDEL;
	`chmod +x $scripts/snp.$chr.sh $scripts/indel.$chr.sh`;
}
