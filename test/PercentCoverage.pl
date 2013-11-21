#!usr/bin/perl

use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
require Subfunctions;

my $ref_file = shift @ARGV;
my $contigs_file = shift @ARGV;
my $gene_name = shift @ARGV;

percentcoverage ($ref_file, $contigs_file, $gene_name);
