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
while (my $line = readline $TARGET_FH) {
	$line =~ /(.*),(.*)/;
	push @targets, $1;
}

print "found targets ".join (",", @targets). "\n";
my $raw_hit_matrix = {};
my $hit_matrix = {};
make_hit_matrix ($blast_file, $raw_hit_matrix);
my @contig_names = ();
my $high_score = process_hit_matrix ($raw_hit_matrix, \@contig_names, \@targets, 300, 200, $hit_matrix);
