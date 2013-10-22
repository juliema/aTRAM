#!usr/bin/perl

$path ='/Users/juliema/Documents/Projects/TRAMiEvoBio/3-EXONSBLAST/*.names.txt.FASTAOUT.fas.fasta';
@array = glob($path);

for $file (@array)
{
    $count++;
#    if ($count < 3)
#    {
  $file =~ s/\S+3-EXONSBLAST\/(\S+)-PA.names.txt.FASTAOUT.fas.fasta/$1/g;
	print "$file\t$count\n";
	system "blastn -query $file\-PA.names.txt.FASTAOUT.fas.fasta -subject BodyLouse.$file\-RA.fasta -out $file.blastn.out";
#    }
}

