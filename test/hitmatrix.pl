#!usr/bin/perl

use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
require Subfunctions;
require SequenceRetrieval;

my $contigs_file = shift @ARGV;
my $target_file = shift @ARGV;
my $out_file = shift @ARGV;

my ($TARGET_FH, $target_fasta) = tempfile(UNLINK => 1);
flattenfasta($target_file, $target_fasta, ",");
my @targets = ();
my @seqs = ();

my $log_fh = *STDOUT;
while (my $line = readline $TARGET_FH) {
	$line =~ />(.*),(.*)/;
	push @targets, $1;
	push @seqs, $2;
}
my $protein = 0;
if (is_protein($seqs[0])) {
	$protein = 1;
}
# make a database from the target so that we can compare contigs to the target.
my (undef, $targetdb) = tempfile(UNLINK => 1);
if ($protein == 1) {
	system_call ("makeblastdb -in $target_file -dbtype prot -out $targetdb.db -input_type fasta", $log_fh);
} else {
	system_call ("makeblastdb -in $target_file -dbtype nucl -out $targetdb.db -input_type fasta", $log_fh);
}

my ($BLAST_FH, $blast_file) = tempfile(UNLINK => 1);

if ($protein == 1) {
	system_call ("blastx -db $targetdb.db -query $contigs_file -out $blast_file -outfmt '6 qseqid sseqid bitscore qstart qend sstart send qlen'", $log_fh);
} else {
	system_call ("tblastx -db $targetdb.db -query $contigs_file -out $blast_file -outfmt '6 qseqid sseqid bitscore qstart qend sstart send qlen'", $log_fh);
}

my $raw_hit_matrix = {};
my $hit_matrix = {};
make_hit_matrix ($blast_file, $raw_hit_matrix);
my @contig_names = ();
my $targetptr = \@targets;

my $high_score = process_hit_matrix ($raw_hit_matrix, \@targets, 100, 200, $hit_matrix, \@contig_names);

print "contig\t" . join ("\t",@targets) . "\ttotal score\n";
print "length\t" . join ("\t",map (length ($_), @seqs)) . "\n";

my @sorted_contigs = sort {$hit_matrix->{$b}->{'total'} <=> $hit_matrix->{$a}->{'total'}} keys $hit_matrix;

if ($out_file) {
	open FH, ">", $out_file;
	my $hashed_seqs = findsequences ($contigs_file, \@sorted_contigs);
	foreach my $contig (@sorted_contigs) {
		print FH ">$contig\n";
		print FH "$hashed_seqs->{$contig}\n";
	}
	close FH;
}

foreach my $contig_name (@sorted_contigs) {
	print "$contig_name\t";
	foreach my $target (@targets) {
		my $score = $hit_matrix->{$contig_name}->{$target};
		if ($score == undef) {
			print "-\t";
		} else {
			print "$score\t";
		}
	}
	my $total = $hit_matrix->{$contig_name}->{"total"};
	print "$total\n";
}

