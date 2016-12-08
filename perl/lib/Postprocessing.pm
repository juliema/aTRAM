#!/usr/bin/env perl
package Postprocessing;
use strict;
use File::Temp qw/ tempfile /;
use System;
use Parsing;
use Configuration;

BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw(percentcoverage);
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw(align_to_ref trim_to_ref align_to_seq);
	Configuration::initialize();
}

sub percentcoverage {
	my $reffile = shift;
	my $contigfile = shift;
	my $outname = shift;
	my $aligner = shift;

	my ($contigs, undef) = align_to_ref ($reffile, $contigfile, $outname, $aligner);

	trim_to_ref ($contigs, "reference");
	return $contigs;
}

sub align_to_ref {
	my $reffile = shift;
	my $contigfile = shift;
	my $outname = shift;
	my $aligner = shift;

	# check that the files we need exist:
	unless ((-e $reffile) && (-e $contigfile)) {
		return undef;
	}

	open REF_FH, "<:crlf", $reffile;
	my $ref = readline REF_FH;
	$ref =~ />(.+)$/;
	my $refname = $1;
	close REF_FH;

	my (undef, $catfile) = tempfile(UNLINK => 1);
	`cat $reffile $contigfile > $catfile`;

	return align_to_seq ($catfile, $refname, $outname, $aligner);
}

sub align_to_seq {
	my $catfile = shift;
	my $refname = shift;
	my $outname = shift;
	my $aligner = shift;

	my $result = 0;
	if ($aligner eq "mafft") {
		$result = run_command (get_bin("mafft"), "$catfile > $outname.align.fasta");
	} else {
		$aligner = "muscle";
		$result = run_command (get_bin("muscle"), "-in $catfile -out $outname.align.fasta");
	}

	if ($result == -1) {
		printlog ("Aligner $aligner could not be found.");
		exit -1;
	}

	if ($result != 0) {
		printlog ("Aligner $aligner failed.");
		return undef;
	}

	my ($seq_hash, undef) = parsefasta ("$outname.align.fasta");
	my $gappedrefseq = delete $seq_hash->{"$refname"};
	$seq_hash->{"reference"} = $gappedrefseq;

	return ($seq_hash);
}

sub trim_to_ref {
	my $seq_hash = shift;
	my $refseq_key = shift;

	my $gappedrefseq = delete $seq_hash->{$refseq_key};

	my $result = 0;
	# as long as there are still gaps in the reference sequence, keep removing the corresponding positions from the contigs.
	while ($gappedrefseq =~ /(\w*)(-+)(.*)/) {
		my $left = $1;
		my $gap = $2;
		my $remainder = $3;
		my $leftlength = length $left;
		my $gaplength = length $gap;

		foreach my $contig (keys %$seq_hash) {
			my $end = $leftlength + $gaplength;
			my $start = $leftlength + 1;
			my ($startseq, $regionseq, $endseq) = split_seq ($seq_hash->{$contig}, $start, $end);
			$seq_hash->{$contig} = "$startseq$endseq";
		}

		$gappedrefseq = "$left$remainder";
	}
	# align the ends of the contig seqs
	foreach my $contig (keys %$seq_hash) {
		if ((length $seq_hash->{$contig}) > (length $gappedrefseq)) {
			# if the contig seq is longer than the gappedrefseq, truncate it
			$seq_hash->{$contig} = substr($seq_hash->{$contig}, 0, length $gappedrefseq);
		} else {
			# if the contig seq is shorter than the gappedrefseq, pad it with gaps.
			$seq_hash->{$contig} .= "-" x ((length $seq_hash->{$contig}) - (length $gappedrefseq));
		}
	}

	$seq_hash->{$refseq_key} = $gappedrefseq;
	return $seq_hash;
}

return 1;
