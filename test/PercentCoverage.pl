#!usr/bin/perl

use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
require Subfunctions;
use File::Temp qw/ tempfile /;

my (undef, $filenames) = tempfile(UNLINK => 1);
system ("ls -l TRAM.*.best.fasta > $filenames");
open NAMES, "<", $filenames;
while (<NAMES>)
{
    if (/(TRAM.(\S+).phum\S+.fasta)/)
	{
		my $file=$1;
		my $gene=$2;
		percentcoverage ("$gene.phum.fasta", "$file", "$gene");
	}
}

sub percentcoverage {
	my $reffile = shift;
	my $contigfile = shift;
	my $gene = shift;

	###### cat files
	my (undef, $catfile) = tempfile(UNLINK => 1);
	system ("cat $reffile $contigfile > $catfile");
	##### muscle alignment
	system ("muscle -in $catfile -out $gene.muscle.fasta");

	my (undef, $fastafile) = tempfile(UNLINK => 1);
	flattenfasta("$gene.muscle.fasta", $fastafile, ",");

	# parse the output file: save the reference as a separate sequence, put the others into an array.
	my $refseq = "";
	my $contigs = {};
	open FH, "<", $fastafile;
	while (my $line = readline FH) {
		$line =~ />(.*?),(.*)$/;
		my $name = $1;
		my $seq = $2;
		if ($name =~ /$gene/) {
			$refseq = $seq;
		} else {
			$contigs->{$name} = $seq;
		}
	}
	close FH1;

	# as long as there are still gaps in the reference sequence, keep removing the corresponding positions from the contigs.
	while ($refseq =~ /(\w*)(-+)(.*)/) {
		my $left = $1;
		my $gap = $2;
		my $remainder = $3;
		my $leftlength = length $left;
		my $gaplength = length $gap;

		foreach my $contig (keys $contigs) {
			$contigs->{$contig} =~ /(.{$leftlength})(.{$gaplength})(.*)/;
			$contigs->{$contig} = "$1$3";
		}

		$refseq = "$left$remainder";
	}

	####### Print out EXON file
	####### Print OUT Table

	open TABLE_FH, ">", "$gene.Table.txt";
	open EXON_FH, ">", "$gene.exons.fasta";

	print EXON_FH ">$gene\n$refseq\n";
	print "contig\ttotal\tpercent\n";
	print TABLE_FH "contig\ttotal\tpercent\n";
	my $total_length = length $refseq;
	foreach my $contig (keys $contigs) {
		print EXON_FH ">$contig\n$contigs->{$contig}\n";
		my $gaps = ($contigs->{$contig} =~ tr/N-//);
		my $total = $total_length - $gaps;
		my $percent = $total / $total_length;
		print "$contig\t$total\t$percent\n";
		print TABLE_FH "$contig\t$total\t$percent\n";
	}

	close TABLE_FH;
	close EXON_FH;
}
