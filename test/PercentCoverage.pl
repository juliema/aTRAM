#!/usr/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
require Subfunctions;

my $ref_file = shift @ARGV;
my $contigs_file = shift @ARGV;
my $gene_name = shift @ARGV;
my $aligner = shift @ARGV;

if ($gene_name eq "") {
	print "Usage: PercentCoverage.pl ref_file contig_file gene_name aligner\n";
	exit;
}

my $contigs = percentcoverage ($ref_file, $contigs_file, $gene_name, $aligner);

my $refseq = delete $contigs->{reference};

open TABLE_FH, ">", "$gene_name.Table.txt";
open EXON_FH, ">", "$gene_name.exons.fasta";

print EXON_FH ">$gene_name\n$refseq\n";
print TABLE_FH "contig\ttotal\tpercent\n";
my $total_length = length $refseq;
foreach my $contig (keys $contigs) {
	print EXON_FH ">$contig\n$contigs->{$contig}\n";
	my $gaps = ($contigs->{$contig} =~ tr/N-//);
	my $total = $total_length - $gaps;
	my $percent = $total / $total_length;
	print TABLE_FH "$contig\t$total\t$percent\n";
}

close TABLE_FH;
close EXON_FH;
