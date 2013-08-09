#!usr/bin/perl

####### There are duplicates, need to get rid of them, it takes time and space

$pathDB='/Users/juliema/Documents/Projects/TRAMiEvoBio/P_schaefBlastDB';


$pathNames= '/Users/juliema/Documents/Projects/TRAMiEvoBio/*.names.txt';

@array = glob("$pathNames");

for $each(@array)
{
    print "$each\n";
    system "blastdbcmd -db $pathDB/s_3_1_sequence_pediculus_blast_database -dbtype nucl -entry_batch $each -outfmt %f -out $each.FASTAOUT.fas";
}     
