#!usr/bin/perl

$seqretrieve=3000;
open TABLE, ">TblastN.table";

open FH, "<Pediculus_hum.peptides.tblastn.max3000.out";
while (<FH>)
{
    if (/Query\S+\s+(\S+)/)
    {
  $count=0;
	$out=$1;
	open OUT, ">$out.names.txt";
    }
    elsif (/(\d+)\s+hits\s+found/)
    {
	$numhits=$1;
	print TABLE "$out\t$numhits\n";
	print  "$out\t$numhits\n";
    }
    elsif (/$out\s+(\S+)/)
    {
	$seq=$1;
	print OUT "$seq\n";
    }
}		


