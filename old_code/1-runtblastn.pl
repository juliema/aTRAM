#!usr/bin/perl
use strict;

my $db_path=shift;
my $queryfile=shift;
my $outfile = $queryfile . ".out";
my $numtargets=3000;

# $path="/Users/juliema/Documents/Projects/TRAMiEvoBio/P_schaefBlastDB";
#
# system "tblastn -db $path/s_3_1_sequence_pediculus_blast_database -outfmt 7 -max_target_seqs $numtargets -query Pediculus_hum.peptides.fas -out tblastn.out";


system "tblastn -db $db_path -outfmt 7 -query $queryfile -out $outfile";

