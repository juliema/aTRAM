#!/usr/bin/env perl
use strict;
use File::Temp qw/ tempfile /;
use FindBin;
use lib "$FindBin::Bin/../lib";
use System;
use Configuration;
use Parsing;

if (@ARGV < 2) {
	print "Usage: blastcontigs.pl contigs_fasta target_fasta [blast_file]\n";
	exit;
}

my $contigs_fasta = shift;
my $target_fasta = shift;
my $blast_file = shift;

Configuration::initialize();

# make a database from the target so that we can compare contigs to the target.
my (undef, $targetdb) = tempfile(UNLINK => 1);
my (undef, $targets) = tempfile(UNLINK => 1);

unless (defined $blast_file) {
	(undef, $blast_file) = tempfile(UNLINK => 1);
}

run_command (get_bin("makeblastdb"), "-in $target_fasta -dbtype nucl -out $targetdb.db -input_type fasta");
run_command (get_bin("tblastx"), "-db $targetdb.db -query $contigs_fasta -out $blast_file -outfmt '6 qseqid sseqid bitscore qstart qend sstart send qlen'");

sortfasta ($target_fasta, $targets, "#");
open FH, "<:crlf", $targets;
my @targets = ();
foreach (my $line = readline FH) {
	$line =~ />(.*?)#/;
	push @targets, $1;
}

my $hit_matrix = {};

make_hit_matrix ($blast_file, $hit_matrix);

print "contig\t" . join ("\t",@targets) . "\ttotal score\n";
foreach my $contig (keys %$hit_matrix) {
	my $contigname = "$contig";
	print "$contigname\t";
	my $total = 0;
	foreach my $target (@targets) {
		if (defined $hit_matrix->{$contig}->{$target}) {
			my $score = abs($hit_matrix->{$contig}->{$target});
			print "$score\t";
			$total += $score;
		} else {
			print "-\t";
		}
	}
	print "$total\n";
}
