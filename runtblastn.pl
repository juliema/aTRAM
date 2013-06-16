#!usr/bin/perl

$numtargets=3000;

$path="/Users/juliema/Documents/Projects/TRAMiEvoBio/P_schaefBlastDB";

system "ls -l *.fas >filenames";

open FH, "<filenames";
while (<FH>)
{
    if (/(\S+).fas/ && ! /Ped_all/)
    {
  $query = $1;
	print "$query\n";
	system "tblastn -db $path/s_3_1_sequence_pediculus_blast_database -outfmt 7 -max_target_seqs $numtargets -query $query.fas -out $query.tblastn.max$numtargets.out ";
    }
}
