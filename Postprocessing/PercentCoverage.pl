#!/usr/bin/env perl

use strict;
use FindBin;
use File::Basename qw(basename);
use lib "$FindBin::Bin/../lib";
use Postprocessing;
use Parsing qw(parsefasta);
use File::Temp qw(tempfile);

my $ref_file = shift @ARGV;
my $contigs_file = shift @ARGV;
my $out_name = shift @ARGV;
my $aligner = shift @ARGV;

if ($out_name eq "") {
	print "Usage: PercentCoverage.pl ref_file contig_file out_name aligner\n";
	exit;
}

if (!(-e $ref_file)) {
	die "Couldn't find reference file $ref_file";
}

if (!(-e $contigs_file)) {
	die "Couldn't find contigs file $contigs_file";
}

my ($reffasta, $reffastaarray) = parsefasta ($ref_file);
my ($fh, $filename) = tempfile();
print $fh ">reference\n";
foreach my $r (@$reffastaarray) {
	print $fh "$reffasta->{$r}\n";
}
close $fh;

my $contigs = percentcoverage ($filename, $contigs_file, $out_name, $aligner);

if (defined $contigs) {
	my $refseq = delete $contigs->{reference};

	my $ref_name = basename($ref_file);
	$ref_name =~ s/\.fa(sta)*//;

	open TABLE_FH, ">", "$out_name.results.txt";
	open EXON_FH, ">", "$out_name.trimmed.fasta";

	print EXON_FH ">$ref_name\n$refseq\n";
	print TABLE_FH "contig\ttotal\tpercent\n";
	my $total_length = length $refseq;
	foreach my $contig (keys %$contigs) {
		print EXON_FH ">$contig\n$contigs->{$contig}\n";
		my $gaps = ($contigs->{$contig} =~ tr/N-//);
		my $total = $total_length - $gaps;
		my $percent = $total / $total_length;
		print TABLE_FH "$contig\t$total\t$percent\n";
	}

	close TABLE_FH;
	close EXON_FH;
}
