#!/usr/bin/perl
use strict;

my $short_read_archive = shift;
my $genes_file = shift;

my $fastq_input = 0;

# if the sra is a fastq file, make it fasta.
if ($short_read_archive =~ /\.f.*q/) { # if it's a fastq file:
	$fastq_input = 1;
	# un-interleave fastq file into fasta:
	open FH, "<", $short_read_archive or die "couldn't open fasta file";

	open OUT_FH, ">", "$short_read_archive.fasta" or die "couldn't create result file";

	my $fs = readline FH;
	while ($fs) {
		$fs =~ s/@/>/;
		print OUT_FH $fs;
		$fs = readline FH;
		print OUT_FH $fs;
		$fs = readline FH;
		# toss quality lines
		$fs = readline FH;
		$fs = readline FH;
	}

	close FH;
	close OUT1_FH;

}

# now we're working with a fasta file for sure.
my $working_sra = $short_read_archive;
if ($fastq_input == 1) {
	$working_sra = "$short_read_archive.fasta";
}

# sort fasta short-read file
system ("bash 3.5-sort_fasta.sh $working_sra");

# un-interleave fasta file into two paired files:
open FH, "<", $short_read_archive or die "couldn't open fasta file";

open OUT1_FH, ">", "$short_read_archive.1.fasta" or die "couldn't create result file";
open OUT2_FH, ">", "$short_read_archive.2.fasta" or die "couldn't create result file";

my $fs = readline FH;
while ($fs) {
	print OUT1_FH $fs;
	$fs = readline FH;
	print OUT1_FH $fs;
	$fs = readline FH;
	print OUT2_FH $fs;
	$fs = readline FH;
	print OUT2_FH $fs;
	$fs = readline FH;
}

close FH;
close OUT1_FH;
close OUT2_FH;

# make the blast db from the first of the paired end files
system ("makeblastdb -in $short_read_archive.1.fasta -dbtype nucl -out $short_read_archive.db");
