#!usr/bin/perl

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../lib";
require Subfunctions;

my $ref_file = 0;
my $contigs_file = 0;
my $output_name = 0;
my $help = 0;

GetOptions ('reference=s' => \$ref_file,
            'input|contigs=s' => \$contigs_file,
            'output=s' => \$output_name,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

open REF_FH, "<", $ref_file;
my $ref = readline REF_FH;
$ref =~ />(.+)$/;
my $refname = $1;
close REF_FH;

my $contigs = percentcoverage ($ref_file, $contigs_file, $refname);
my $refseq = delete $contigs->{reference};

open EXON_FH, ">", "$output_name.exons.fasta";

print EXON_FH ">$refname\n$refseq\n";
foreach my $contig (keys $contigs) {
	print EXON_FH ">$contig\n$contigs->{$contig}\n";
}

close EXON_FH;
