# THIS IS OBSOLETED: use module instead

# given two fasta files corresponding to paired ends of short reads, find the specified sequences and then output them into a single output file.
#!/usr/bin/perl
use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
require pairedsequenceretrieval;

if (@ARGV < 3) {
	die "Usage: 2.5-pairedsequenceretrieval.pl shortreads.#.fasta sequencelist outfile\n";
}

my $fastafile = shift;
my $sequencelist = shift;
my $outfile = shift;

pairedsequenceretrieval ($fastafile, $sequencelist, $outfile);
