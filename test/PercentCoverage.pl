#!usr/bin/perl

use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
require Subfunctions;
use File::Temp qw/ tempfile /;

my (undef, $filenames) = tempfile(UNLINK => 1);
system ("ls -l TRAM.*.best.fasta > $filenames");
open NAMES, "<", $filenames;
while (<NAMES>)
{
    if (/(TRAM.(\S+).phum\S+.fasta)/)
	{
		my $file=$1;
		my $gene=$2;
		percentcoverage ("$gene.phum.fasta", "$file", "$gene");
		system ("cat $gene.Table.txt");
	}

}

