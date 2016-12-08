#!/usr/bin/env perl
use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Sequenceretrieval;

if (@ARGV < 3) {
    die "Usage: findsequences.pl fastafile sequencelist outfile\n";
}

my $fastafile = shift;
my $sequencelist = shift;
my $outfile = shift;

unless (-e $sequencelist) {
	die "File $sequencelist does not exist.\n";
}

open FH, "<:crlf", $sequencelist;
my @sequences = ();
while (my $line = readline FH) {
	chomp $line;
	push @sequences, $line;
}
close FH;

my $results = findsequences ($fastafile, \@sequences);

open OUT_FH, ">", $outfile;

for (my $i=0; $i<@sequences; $i++) {
	print OUT_FH ">$sequences[$i]\n$results->{$sequences[$i]}\n";
}

close OUT_FH;
