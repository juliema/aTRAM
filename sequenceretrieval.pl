#!usr/bin/perl                                                                                                                                                

system "ls -l *blastn3.out >filenames";

open FH, "<filenames";
while (<FH>)
{
    if (/(\S+).blastn3.out/)
    {
        %hash=();
        $count++;
        $gene=$1;
        print "$gene\t$count\n";
        open OUT, ">$gene.names.3.out";
        open FH1, "<$gene.blastn3.out";
        while (<FH1>)
        {
            if (/Query\S+\s+(\S+)/)
            {
#               $count=0;                                                                                                                                     
                $out=$1;
            }
            elsif (/^$out\s+(\S+)/)
            {
                $seq=$1;
                if (! exists $hash{$seq})
                {
                    print OUT "$seq\n";
                }
            }
        }
    }
}

close FH;
close FH1;

$pathDB='/Users/juliema/Documents/1-Projects/P_schaefBlastDB';

$pathNames= '/Users/juliema/Documents/1-Projects/TRAMiEvoBio/1THOUSANDGENES/*.names.3.out';

@array = glob($pathNames);

for $each(@array)
{
    $count++;
    print "$each\t$count\n";
    system "blastdbcmd -db $pathDB/s_3_1_sequence_pediculus_blast_database -dbtype nucl -entry_batch $each -outfmt %f -out $each.FASTAOUT.3.fas";
}

