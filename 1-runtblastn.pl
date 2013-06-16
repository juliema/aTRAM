#!usr/bin/perl

$numtargets=3000;

### I donnt think the numtargests worked 


$path="/Users/juliema/Documents/Projects/TRAMiEvoBio/P_schaefBlastDB";

system "tblastn -db $path/s_3_1_sequence_pediculus_blast_database -outfmt 7 -max_target_seqs $numtargets -query Pediculus_hum.peptides.fas -out tblastn.out";

