#!usr/bin/perl


######  remove duplicate files
######  done

$pathfiles = '/Users/juliema/Documents/1-Projects/TRAMiEvoBio/1THOUSANDGENES/*.FASTAOUT.3.fas';
$pairedend = '/Users/juliema/Documents/Pediculus_schaeffi_gDNA/pediculus/s_3_2_sequence_pediculus.txt'; 



@array=glob("$pathfiles");

for $file (@array)
{
    $flag=0;
    $numseqs=0;
    $seqsfound=0;
    $countfiles++;
    $unique=0;
    %seqhash=();
    %foundseqs=();
    @seqarray=();
    $countnotfound=0;
    print "$countfiles\t$file\n";
    $file =~ s/(\S+).FASTAOUT.3.fas/$1/g;
    open FH, "<$file.FASTAOUT.3.fas";
    open OUT, ">$file.BOTH.3.fas";
    while (<FH>)
    {
	print OUT;
	if (/^>lcl\|(.*)\/1/)
	{
	    $numseqs++;
	    $seq=$1;
	    push @seqarray, $seq;
	    if (! exists $seqhash{$seq})
	    {
		$unique++;
		$seqhash{$seq}=1;
	    }
	}
	    
    }
    close FH;
    print "the number of sequences are $numseqs\n";
    if ($numseqs > 1) 
    {
	open FH1, "<$pairedend";
	while (<FH1>)
	{
	    if (/^@(.*)\/2/)
	    {
		$name=$1;
#		print "$name\n";
		if (exists $seqhash{$name})
		{
		    print OUT ">$name/2\n";
		    $flag=1;
		    $foundseqs{$name}=1;
		    $seqsfound++;
#		    print "found seq $seqsfound\n";
		}
	    }
	    elsif ($flag == 1)
	    {
		$line=$_;
		print OUT "$line";
		$flag = 0;
		if ($seqsfound == $numseqs) { close FH1; }
	    }
	}
    }
    for $seq (@seqarray)
    {
	if (! exists $foundseqs{$seq})
	{
	    $new++;
	    print "This one was not found $seq\n";
	    $countnotfound++;
	}
    }
    print "seqs $numseqs\t\tarray seqs $new\tunique $unique\tseqsfound\t$seqsfound\tseqs not found\t$countnotfound\n";
}
	    
