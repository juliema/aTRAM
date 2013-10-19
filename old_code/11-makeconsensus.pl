#!usr/bin/perl


## 7.29.13
## need files called ambiguity.txt
## script Might be reading Ns properly now:  need to check
## now all ouTput files are the same length !!! BUT NOT TAKING INTO ACCOUNT GAPS IN BODY LOUSE SEQUENCE### ALIGNMENTS ARE ###  OFF

open AMBIG, "<amiguity.txt";
while (<AMBIG>)
{
    if (/(\S+)\s+(\S+)/ && ! /^\#/)
    {
	$ambighash{$1}=$2;
#	print "$1\t$2\n";
    }
}


$path ='/Users/juliema/Documents/1-Projects/TRAMiEvoBio/6-CONSENSUS/*.Chimp.fasta';
#$path ='/Users/juliema/Documents/Projects/TRAMiEvoBio/6-CONSENSUS/FILESTOCHECK/*.Chimp.fasta';

@filearray = glob($path);
for $file (@filearray)
{
    $newlength=0;
    $diff = 0;
#    print "$file\n\n";
#    $file =~ s/\S+FILESTOCHECK\/(\S+).Chimp.fasta/$1/g;
    $file =~ s/\S+CONSENSUS\/(\S+).Chimp.fasta/$1/g;
    $countfiles++;
    print "$file\t$countfiles\n";
    open FH, "<$file.Chimp.fasta";
    open OUT, ">$file.Chimp.conseq.fasta";
    print OUT ">$file.chimp.consensus\n";
    @bigarray=();
    $pos=0;
    while (<FH>)
    {
	if (/\S+/ && ! /^>/)
	{
	    $seq=$_;
	    chomp $seq;
	    $length=length($seq);
	    if ($length > $newlength)
	    {
		$newlength = $length;
	    }
	}
    }
    close FH;
    open FH, "<$file.Chimp.fasta";
    while (<FH>)
    {
	if (/\S+/ && ! /^>/)
	{
	    $seq=$_;
	    chomp $seq;
	    $length=length($seq);
	    @array = split(//, $seq);
	    for $nuc (@array)
	    {
		push @{$bigarray[$pos]}, $nuc;
	    }
	    if ($length < $newlength)
	    {
		$diff = $newlength - $length;
		print "$file\tThe difference is $diff\n";
		for (1..$diff)
		{
		    $nuc = '-';
		    push @{$bigarray[$pos]}, $nuc;
		}
	    } 
	    $pos++;
	}
    }
    for (0..$newlength-1)
    {
	%nuchash=();
	$numitems = 0;
	$nucpos=$_;
	@itemsarray=();
	for (0..$pos)
	{
	    $seqnum=$_;
	    ##### ask if doesnt exist in big array move on.
	    $nuc = $bigarray[$seqnum][$nucpos];
#	    print "$nuc\t";
	    if (! exists $nuchash{$nuc})
	    {
		$numitems++;
		$nuchash{$nuc}=1;
		push @itemsarray, $nuc;
	    }
	    elsif (exists $nuchash{$nuc})
	    {
		$nuchash{$nuc}++;
	    }
#	    print "$nuchash{$nuc}\n";
	}
	$numitems = $numitems-1;
##    print "The number of items is $numitems\n";
##    for $each (@itemsarray) {print "$each\t";} print "\n";
	if ($numitems == 1)
	{
	    $nuc = $itemsarray[0];
	    if ($nuc =~ m/-/)
	    {
		print OUT "N";
#		print "N\t";
	    }
	    else   
	    { 
		print OUT "$nuc";
#		print "$nuc\t";
	    } 
	}
	elsif ($numitems >  1) 
	{
	    if ($numitems == 2)
	    {
		$newlength=0;
		$flag=0;
		for $each (@itemsarray)
		{
		    if ($each =~ m/-/)
		    {
			$flag=1;
		    }
		    elsif ( $each =~ m/N/)
		    {
			$flag=1;
		    }
		}
		if ($flag == 1) 
		{
		    for $each(@itemsarray) 
		    { 
			if ($each =~  m/-/) {} 
#			elsif ($each =~  m/N/) {} 
			else 
			{ 
			    print OUT "$each";  
#                           print "$each\t"; 
			}
		    }
		}
		if ($flag == 0) 
		{    
		    $set=();
		    $ambig=();
		    for $each(@itemsarray)
		    { 
			$length = $nuchash{$each}; 
			if ($length > $newlength) 
			{
			    $ambig=0;
			    $newlength = $length; 
			    $nuctoprint = $each;
			    $set=$each;
			}
			elsif ($length == $newlength)
			{
			    $ambig=1;
##			print "length $each $set\n";
			    $set = $set.$each;
##			print "$set\n";
			}
		    }
		    if ($ambig == 1)
		    {
			$set=uc($set);
			print OUT "$ambighash{$set}";
#			print  "$ambighash{$set}\t";
		    }
		    elsif ($ambig == 0)
		    {
			print OUT "$nuctoprint";
#			print "$nuctoprint\t";
		    }
		}
	    }
	    elsif ($numitems > 2 )
	    {
		$set=();
		$newlength=0;
		for $each(@itemsarray)
		{
		    if ($each =~ m/-/) {}
		    elsif ($each =~ m/N/) {}
		    else
		    {
#			print "$each\t";
			$length=$nuchash{$each};
#			print "\n$nucpos\t$each\t$length\n";
			if ($length > $newlength)
			{
			    $ambig=0;
			    $newlength = $length; 
			    $nuctoprint = $each;
			    $set=$each;
#			    print "nuc to print = $nuctoprint position $nucpos\n";
			}
			elsif ($length == $newlength)
			{
##			print "ambig length = new length $each\n";
			    $ambig=1;
##			print "length $set $each\n";
			    $set = $set . $each;
#			print "$nucpos\tAmbigs are $set\n";
			}
		    }
		}
		if ($ambig == 1)
		{
		    $set=uc($set);
		    $amlen = length($set);
		    if ($amlen == 1)  
		    {
			print OUT "$set";  
#			print "$nucpos\t$set\n";
		    }
		    else
		    {
			print OUT "$ambighash{$set}";
#			print "$ambighash{$set}\t";
		    }    
		}
		elsif ($ambig == 0)
		{
#		    print  "$nuctoprint\t";
		    print OUT "$nuctoprint";
		}
	    }
	}
    }
#    print "\n";
}
