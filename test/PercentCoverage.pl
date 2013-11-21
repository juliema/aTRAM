#!usr/bin/perl

use strict;
use warnings;

my $name;
my $seq;
my $totalseq;
my $flag;
my @location;
my $humlength;
my $nuc;
my @array;
my $total;
my $position;
my $percent;
my %percenthash;
my %totalhash;
my @namearray;
my $gene;
my %humanhash;
my $totalref;
my $file;

system "ls -l TRAM.*.Best.fasta >filenames";
open NAMES, "<filenames";
while (<NAMES>)
{
    if (/(TRAM.(\S+).phum\S+.fasta)/)
    {
	$file=$1;
	$gene=$2;
	open TABLE, ">$gene.Table.txt";
	open EXON, ">$gene.exons.fasta";
###### cat files
	system "cat $gene.phum.fasta $file >$gene.TRAMOUT.fasta";
##### muscle alignment
	system "./muscle3.8.31_i86darwin64 -in $gene.TRAMOUT.fasta -out $gene.TRAMOUT.muscle.fasta";
#### edit FASTA 
	open FH, "<$gene.TRAMOUT.fasta";
	open ED, ">$gene.TRAMOUT.ed.fasta";
	while (<FH>)
	{
	    if (/^>/)
	    {
		$name=$_;
		print ED "$totalseq\n$name";
		$totalseq=();
	    }
	    elsif (/(\S+)/)
	    {
		$seq=$1;
		$totalseq=$totalseq . $seq;
	    }
	}
	print ED "$totalseq\n";
	close ED;
	close FH;
######################################
#### calculate Percentage ###########
#####################################

### Get Length of Reference 
############################
	open FH1, "<$gene.TRAMOUT.ed.fasta";
	while (<FH1>)
	{
	    if ($flag == 1)
	    {
		$totalref=();
		@location =();
		$humlength=0;
		$seq=$_;
		$flag++;
		@array = split(//,$seq);
		for $nuc (@array)
		{
		    if ($nuc ne '-')
		    {
			push @location, ($humlength);
			$totalref++;
		    }
		    $humlength++;
		}
		print  "$gene human $totalref\n";
	    }
	    if (/phum/)
	    {
		$name=$_;
		chomp $name;
		$flag=1;
	    }
	}
	$humanhash{$gene}=$totalref;
	close FH1;
####### Calculate lengths of all sequences 
####### Print out EXON file
####### Print OUT Table
###########################################
	open FH2, "<$gene.TRAMOUT.ed.fasta";
	while (<FH2>)
	{
	    if (/^>(\S+)/)
	    {
		$name=$1;
		print EXON ">$name\n";
		push @namearray, $name;
		print "The name is $name\n";
	    }
	    elsif (/(\S+)/)
	    {
		$total=();
		$seq=$1;
		@array = split(//,$seq);
		for $position (@location)
		{
		    $nuc = $array[$position];
		    print EXON "$nuc";
		    if ($nuc ne "N")
		    {
			if ($nuc ne "-")
			{
			    $total++;
			}
		    }
		}
		print EXON "\n";
	    }
	    $percent = $total/$totalref;
	    $percenthash{$name}=$percent;
	    $totalhash{$name}=$total;
	}
	for $name(@namearray)
	{
	    print  "$name\t$totalhash{$name}\t$percenthash{$name}\n";
	    print TABLE "$name\t$totalhash{$name}\t$percenthash{$name}\n";
	}		
    }
}    





