#!usr/bin/perl

$path = '/Users/juliema/Documents/Projects/TRAMiEvoBio/5-MUSCLE/*blastn.out.fasta';

@array = glob($path);
for $file (@array)
{
    $count++;
    $fullseq = ();
    $file =~ s/\S+\/5-MUSCLE\/(\S+).blastn.out.fasta/$1/g;
#    print "$file\n";
    print "$count\n";
    open FH, "<BodyLouse.$file\-RA.fasta";
    while (<FH>)
    {
  if (/^>body\_louse\_(\S+)-RA/)
	{
	    $name=$1;
	}
        else
        {
	    $seq=$_;
	    chomp $seq;
	 }
        $fullseq=$fullseq . $seq;
     }
     open FH1, "<$file.blastn.out.fasta";
     while (<FH1>)
     {
	 if (/^>(\S+)/)
         {
	     $contig=$1;
	     open OUT, ">$file.$contig.toalign.fas";
	     print OUT ">BodyLouse.$name\n$fullseq\n";
	     print OUT ">$contig\n";
	 }
	 else
         {
	     print OUT;
	 }
     }
}
