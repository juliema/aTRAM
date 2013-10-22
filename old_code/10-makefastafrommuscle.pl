#!usr/bin/perl

$path ='/Users/juliema/Documents/Projects/TRAMiEvoBio/6-CONSENSUS/*.muscle.out';

@array=glob($path);
for $file (@array)
{
#    $count++;
#    if ($count < 100)
#    {
  $file =~ s/\S+CONSENSUS\/(\S+)\.\d+_\d+.muscle.out/$1/g;
	if (! exists $genehash{$file})
	{
	    $countfiles++;
	    print "$file\t$countfiles\n";
	    $genehash{$file}=1;
	    open OUT, ">$file.Chimp.fasta";
	    system "ls -l $file*.muscle.out >filenames";
	    open FH, "<filenames";
	    while (<FH>)
	    {
		if (/$file.(\S+).muscle.out/)
		{
		    $read=$1;
		    $fullseq=();
		    open FH1, "<$file.$read.muscle.out";
		    while (<FH1>)
		    {
			if (/^>/ && ! /BodyLouse/)
			{
			    $flag=1;
			    print OUT;
			}
			elsif (/^>BodyLouse/)
			{
			    $flag =0;
			}
			elsif ($flag == 1)
			{
			    $seq=$_;
			    chomp $seq;
			    $fullseq=$fullseq.$seq;
			}
		    }
		    print OUT "$fullseq\n";
		}
	    }
	}
#    }
}

