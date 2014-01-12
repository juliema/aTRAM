#!/usr/bin/perl
use strict;
use File::Temp qw/ tempfile /;
use FindBin;
use lib "$FindBin::Bin/../lib";
require Subfunctions;

my $contigs_fasta = shift;
my $target_fasta = shift;
my $blast_file = shift;

if ($gene_name eq "") {
	print "Usage: blastcontigs.pl contigs_fasta target_fasta blast_file\n";
	exit;
}


my $blast_asn = "$blast_file.asn";

# make a database from the target so that we can compare contigs to the target.
my (undef, $targetdb) = tempfile(UNLINK => 1);
my (undef, $targets) = tempfile(UNLINK => 1);

unless ($blast_file) {
	(undef, $blast_file) = tempfile(UNLINK => 1);
	(undef, $blast_asn) = tempfile(UNLINK => 1);
}

system_call ("makeblastdb -in $target_fasta -dbtype nucl -out $targetdb.db -input_type fasta");
system_call ("tblastx -db $targetdb.db -query $contigs_fasta -out $blast_asn -outfmt 11");

system_call ("blast_formatter -archive $blast_asn -out $blast_file -outfmt '6 qseqid sseqid bitscore qstart qend sstart send qlen'");


sortfasta ($target_fasta, $targets, "#");
open FH, "<", $targets;
my @targets = ();
foreach (my $line = readline FH) {
	$line =~ />(.*?)#/;
	push @targets, $1;
}

my $hit_matrix = {};

make_hit_matrix ($blast_file, $hit_matrix);

print "contig\t" . join ("\t",@targets) . "\ttotal score\n";
foreach my $contig (keys $hit_matrix) {
	my $contigname = "$contig";
	print "$contigname\t";
	my $total = 0;
	foreach my $target (@targets) {
		if ($hit_matrix->{$contig}->{$target} == undef) {
			print "-\t";
		} else {
			my $score = abs($hit_matrix->{$contig}->{$target});
			print "$score\t";
			$total += $score;
		}
	}
	print "$total\n";
}
