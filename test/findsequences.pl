#!/usr/bin/perl
use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
require Sequenceretrieval;

if (@ARGV < 3) {
    die "Usage: 6.5-findsequences.pl fastafile sequencelist outfile\n";
}

my $fastafile = shift;
my $sequencelist = shift;
my $outfile = shift;

findsequences ($fastafile, $sequencelist, $outfile);
