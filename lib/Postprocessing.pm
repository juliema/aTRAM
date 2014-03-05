#!/usr/bin/env perl
package Postprocessing;
use strict;
use File::Temp qw/ tempfile /;
use System;
use Parsing;

BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw(percentcoverage);
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw();
}

sub percentcoverage {
	my $reffile = shift;
	my $contigfile = shift;
	my $outname = shift;
	my $aligner = shift;

	# check that the files we need exist:
	unless ((-e $reffile) && (-e $contigfile)) {
		return undef;
	}

	my (undef, $catfile) = tempfile(UNLINK => 1);
	my $result = 0;
	$result = system_call ("cat $reffile $contigfile > $catfile");

	if ($aligner eq "mafft") {
		$result = system_call ("mafft $catfile > $outname.align.fasta");
	} else {
		$aligner = "muscle";
		$result = system_call ("muscle -in $catfile -out $outname.align.fasta");
	}

	if ($result == -1) {
		printlog ("Aligner $aligner could not be found.");
		exit -1;
	}

	if ($result != 0) {
		printlog ("Aligner $aligner failed.");
		return undef;
	}

	open REF_FH, "<", $reffile;
	my $ref = readline REF_FH;
	$ref =~ />(.+)$/;
	my $refname = $1;
	close REF_FH;

	# parse the output file: save the reference as a separate sequence, put the others into an array.
	my ($contigs, $taxa) = parsefasta ("$outname.align.fasta");
	my $refseq = delete $contigs->{"$refname"};

	# as long as there are still gaps in the reference sequence, keep removing the corresponding positions from the contigs.
	while ($refseq =~ /(\w*)(-+)(.*)/) {
		my $left = $1;
		my $gap = $2;
		my $remainder = $3;
		my $leftlength = length $left;
		my $gaplength = length $gap;

		foreach my $contig (keys %$contigs) {
			my $end = $leftlength + $gaplength;
			my $start = $leftlength + 1;
			my ($startseq, $regionseq, $endseq) = split_seq ($contigs->{$contig}, $start, $end);
			$contigs->{$contig} = "$startseq$endseq";
		}

		$refseq = "$left$remainder";
	}
	# align the ends of the contig seqs
	foreach my $contig (keys %$contigs) {
		if ((length $contigs->{$contig}) > (length $refseq)) {
			# if the contig seq is longer than the refseq, truncate it
			$contigs->{$contig} = substr($contigs ->{$contig}, 0, length $refseq);
		} else {
			# if the contig seq is shorter than the refseq, pad it with gaps.
			$contigs->{$contig} .= "-" x ((length $contigs->{$contig}) - (length $refseq));
		}
	}

	$contigs->{"reference"} = $refseq;
	return $contigs;
}

return 1;
