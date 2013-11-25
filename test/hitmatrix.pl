#!usr/bin/perl

use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
require Subfunctions;

my $blast_file = shift @ARGV;
my $target_file = shift @ARGV;

my ($TARGET_FH, $target_fasta) = tempfile(UNLINK => 1);
flattenfasta($target_file, $target_fasta, ",");
my @targets = ();
my @seqs = ();
while (my $line = readline $TARGET_FH) {
	$line =~ />(.*),(.*)/;
	push @targets, $1;
	push @seqs, $2;
}

my $raw_hit_matrix = {};
my $hit_matrix = {};
make_hit_matrix ($blast_file, $raw_hit_matrix);
my @contig_names = ();
my $targetptr = \@targets;

my $high_score = process_hit_matrix ($raw_hit_matrix, \@targets, 100, 200, $hit_matrix, \@contig_names);

print "contig\t" . join ("\t",@targets) . "\ttotal score\n";
print "length\t" . join ("\t",map (length ($_), @seqs)) . "\n";

foreach my $contig_name (@contig_names) {
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

