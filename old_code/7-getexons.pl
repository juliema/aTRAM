#!usr/bin/perl

$path = '/Users/juliema/Documents/Projects/TRAMiEvoBio/4-ALIGNMENTS/*.blastn.out';

@array = glob($path);
for $file (@array)
{

    $count++;
#    if ($count < 3)
#    {
  print "$count\n";
	$file =~ s/(\S+).blastn.out/$1/g;
	open OUT, ">$file.blastn.out.fasta";
	open FH, "<$file.blastn.out";
	$nohits=0;
	while (<FH>)
	{
	    if (/Query=\s+(\S+)/)
	    {
		$line++;
		$flag = 1;
		$query=$1;
		$fullseq=();
		$first = 1;
	    }
	    if ($flag == 1)
	    {	
		if (/No\shits\sfound/)
		{	    
		    $nohits++;
		    $flag = 0;
		    $go=0;
		    $first=0;
		}
		elsif (/Strand=Plus\/(\S+)/)
		{
		    $orientation=$1;
		    print "$query\t$orientation\n";
		}
		elsif (/Query\s+\d+\s+(\S+)/)
		{
		    $seq=$1;
		    $fullseq = $fullseq . $seq;
		    $go=1;
		}
		if ($first == 1)
		{
		    if (/Sbjct\s+(\d+)/)
		    {
			$start=$1;
			$first=0;
			print "Start\t$start\n";
		    }
		}
	    }
	    if (/Effective\ssearch\s/)
	    {
		if ($flag == 1)
		{
		    $flag = 0;
		    $first = 1;
		    $counthits++;
		}
		if ($go == 1)
		{
		    if ($orientation eq 'Plus')
		    {
			print OUT ">$query\_$start\n$fullseq\n";
		    }
		    elsif ($orientation eq 'Minus')
		    {
			$rev=reverse($fullseq);
			$rev =~ tr/ACGT/TGCA/;
			print OUT ">$query\_$start\n$rev\n";
		    }
		}
	    }
	}
	print "There were $counthits matches and $nohits non matches\n";
#    }
}
