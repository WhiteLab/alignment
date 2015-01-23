#!/usr/bin/perl
use strict;
use warnings;

my %files = archiveLookup();

#print "BionimbusID\tDate\tMachine\tRun\tBarCode\tLane\ttota_reads\tpost_align_readst\tfraction_aligned\tRmDupAlignedReads\tFracDup\tOnTargetBases\tReadLength\tFracOnTarget\tTargetEnrichment\tBasesCov>=8x\tTargetBases\t%Cov>=8x\n";
print "BionimbusID\tDate\tMachine\tRun\tBarCode\tLane\tread_length\ttotal_reads\tpost_align_reads\tfraction_aligned\tpost_rmdup_reads\tfraction_rmduped\ttarget_size\taligned_bp_ot\tfraction_aligned_bp_ot\tfraction_sequenced_bp_ot\taverage_ot_cov\tcoverage1\tcoverage2\tcoverage8\tfraction_coverage1\tfraction_coverage2\tfraction_coverage8\taligned_target_enrichment\tsequenced_target_enrichment\n";
my %destinations;
my %runs;
my $GENOME_SIZE = 3137161264;

my @samples = sort keys %files;
foreach my $sample (@samples) {
	print STDERR "$sample\n";
	my @set = sort keys %{$files{$sample}};
    
}
my %data = directory_search(".");
output_stats(\%data);

sub output_stats
{
    my $data_ref = $_[0];
    my %data = %$data_ref;
    my @samples = sort (keys %data);
    foreach my $sample (@samples) {

		my @annot = split (/\_/, $sample);
		my $db_id = $annot[0];
	    my $date = $annot[1];
	    my $machine = $annot[2];
	    my $run = $annot[3];
	    my $barcode = $annot[4];
	    my $lane = $annot[5];

	    unless ($data{$sample}{LEN}) {
			$data{$sample}{LEN} = 100;
	    }
	    unless ($data{$sample}{READS}) {
			$data{$sample}{READS}="NA";
	    }
	    if ($data{$sample}{RAW}{MAP}) {
		    # Print DB ID, Run date, machine SN, run #, flow cell barcode, lane
			print "$db_id\t$date\t$machine\t$run\t$barcode\t$lane\t";

			# get Total reads, Aligned reads, rmdup aligned reads, Fraction duplicate reads
			my $ali_frac = sprintf ("%.4f", $data{$sample}{RAW}{MAP}/$data{$sample}{RAW}{TOTAL});
			my $dup_frac = sprintf ("%.4f", 1 - ($data{$sample}{RMDUP}{MAP} / $data{$sample}{RAW}{TOTAL}));

			# get # on target bases, Read length, Fraction on target bases, On target Enrichment
			my $total_bp = $data{$sample}{RAW}{TOTAL} * $data{$sample}{LEN};
			my $aligned_rmdup_bp = $data{$sample}{RMDUP}{MAP} * $data{$sample}{LEN};
			my $on_target_reads = $data{$sample}{COV}{TOTAL} / $data{$sample}{LEN};
			my $on_target_bp = $data{$sample}{COV}{TOTAL};
			my $frac_aligned_bp_ot = sprintf ("%.4f", $on_target_bp / $aligned_rmdup_bp);
			my $frac_sequenced_bp_ot = sprintf ("%.4f", $on_target_bp / $total_bp);
			my $average_ot_coverage = sprintf ("%.2f", $on_target_bp / $data{$sample}{COV}{TARGET});
			my $on_target_frac = sprintf ("%.4f", $on_target_reads / $data{$sample}{RMDUP}{MAP} );
			
			# calculate enrichment
			
			my $frac_sizes = sprintf ("%.4f", $GENOME_SIZE / $data{$sample}{COV}{TARGET});
			my $frac_aligned_reads = sprintf ("%.4f", $aligned_rmdup_bp / $on_target_bp);
			my $aligned_target_enrichment = sprintf ("%.2f", $frac_sizes / $frac_aligned_reads);

			my $frac_sequenced_reads = sprintf ("%.4f", $total_bp / $on_target_bp);
			my $sequenced_target_enrichment = sprintf ("%.2f", $frac_sizes / $frac_sequenced_reads);

#HEADER: read_length total_reads post_align_reads fraction_aligned post_rmdup_reads fraction_rmduped target_size aligned_bp_ot fraction_aligned_bp_ot fraction_sequenced_bp_ot average_ot_cov coverage1 coverage2 coverage8 fraction_coverage1 fraction_coverage2 fraction_coverage8 aligned_target_enrichment sequenced_target_enrichment\n";
			print "$data{$sample}{LEN}\t";			# read_length
			print "$data{$sample}{RAW}{TOTAL}\t";	# total_reads
			print "$data{$sample}{RAW}{MAP}\t";		# post_align_reads
			print "$ali_frac\t";					# fraction_aligned
			print "$data{$sample}{RMDUP}{MAP}\t";	# post_rmdup_reads
			print "$dup_frac\t";					# fraction_rmduped
			print "$data{$sample}{COV}{TARGET}\t";	# target_size
			print "$data{$sample}{COV}{TOTAL}\t";	# aligned_bp_ot
			print "$frac_aligned_bp_ot\t";			# fraction_aligned_bp_ot
			print "$frac_sequenced_bp_ot\t";		# fraction_sequenced_bp_ot
			print "$average_ot_coverage\t";			# average_ot_cov
			for('1','2','8') {            			# cov 1x 2x 8x
				print $data{$sample}{COV}{"RS$_"} . "\t";
			}
			for('1','2','8') {						# frac cov 1x 2x 8x
				print sprintf ("%.4f", $data{$sample}{COV}{"RS$_"} / $data{$sample}{COV}{TARGET}) . "\t";
			}
			print "$aligned_target_enrichment\t";
			print "$sequenced_target_enrichment\n";
	    }
    }
}
sub directory_search
{
    my %data;
    my @dir = `ls $_[0]`;
    foreach my $file (@dir) {
		chomp $file;
		if ($file =~ /(.+[1-8])\.rmdup.srt.bam.flagstats$/ || $file =~ /(.+[1-8]).rmdup.accepted_hits.bam.flagstats$/) {
			my $sample = $1;
			%{$data{$sample}{RMDUP}} = parseFS("$_[0]/$file");
			my @name = split (/_/, $sample);
			$data{$name[0]}{RUNS}{$_[0]}{$name[5]}++;
				
		}
		elsif ($file =~ /(.+[1-8])\.srt\.bam\.flagstats$/ || $file =~ /(.+[1-8]).accepted_hits.bam.flagstats$/) {
			my $sample = $1;
			%{$data{$sample}{RAW}} = parseFS("$_[0]/$file");
		}
		elsif ($file =~ /(.+[1-8]).(genome|exome|capture)\.hist$/ || $file =~ /(.+[1-8]).rnaSeq.covWhist$/) {
			my $sample = $1;
			($data{$sample}{COV}{TOTAL}, $data{$sample}{COV}{RS1}, $data{$sample}{COV}{RS2}, $data{$sample}{COV}{RS8}, $data{$sample}{COV}{TARGET}) = parseCoverage("$_[0]/$file");
		}
		elsif ($file =~ /(.+[1-8])_[12]_sequence.txt.gz$/ || $file =~ /(.+[1-8])_sequence.txt.gz$/) {
			my $sample = $1;
			my $temp = `gzip -dc $_[0]/$file | head -n 2 | tail -n 1`;
			chomp $temp;
			my $length = length($temp);
			$data{$sample}{LEN} = $length;
		}
		
    }
    return(%data);
}
sub parseFS
{
    my %data;
    my ($fs)  = @_;
    open IN, $fs or die "Can't open $fs!\n";
	my $line = <IN>;
	chomp $line;
	if($line =~ /^(\d+)\s/) {
		$data{TOTAL} = $1;
	}
	$line = <IN>;
	$line = <IN>;
	chomp $line;
	if($line =~ /^(\d+).*\((\d+\.\d+)/) {
		$data{MAP} = $1;
		$data{PCT} = $2;
	}
	return(%data);
}
sub parseCoverage
{
    my $total_bases = 0;
    my $rs1_bases = my $rs2_bases = my $rs8_bases = 0;
    my $target_size = 0;
    my ($fs)  = @_;
    open IN, $fs or die "Can't open $fs!\n";
	my $target_size_flag = 0;
    while (<IN>) {
		chomp $_;
		my @temp = split (/\s+/, $_);
		if(!$target_size_flag){
			$target_size = $temp[3];
			$target_size_flag = 1;
		}
		$total_bases = $total_bases + $temp[1]*$temp[2];
		$target_size = $temp[3];
		if ($temp[1] >= 1) {
			$rs1_bases = $rs1_bases + $temp[2];
		}
		if ($temp[1] >= 2) {
			$rs2_bases = $rs2_bases + $temp[2];
		}
		if ($temp[1] >= 8) {
			$rs8_bases = $rs8_bases + $temp[2];
		}
    }
    return($total_bases, $rs1_bases, $rs2_bases, $rs8_bases, $target_size);
}
sub annotIn
{
    my %annotation;
    open IN, $_[0] or die "Can't open $_[0]!\n";
    while (<IN>) {
	chomp $_;
	my @temp = split (/\t/, $_);
	my $db_id = $temp[0];
	$annotation{$db_id}{EXP_TYPE} = $temp[3];
	$annotation{$db_id}{SAMPLE_NAME} = $temp[1];
	$annotation{$db_id}{USER} = $temp[6];
	$annotation{$db_id}{PROJECT} = $temp[2];
    }
    return(%annotation);
}
sub archiveLookup
{
    my %files;
	my @dir = `ls .`;
	foreach my $file (@dir) {
		chomp $file;
		if ($file =~ /(.+)_sequence.txt(.gz)*/) {
			my @temp = split (/_/, $file);
			my $db_id = $temp[0];
			my $date = $temp[1];
			my $mid = $temp[2];
			my ($year, $month, $day) = (0,0,0);
			if ($date =~ /(\d\d)(\d\d)(\d\d)/) {
				$year = $1;
				$month = $2;
				$day = $3;
			}
			$files{$db_id}{$file}++;
		}
	}
    return(%files);
}

