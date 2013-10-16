#!usr/bin/perl                                                                                                                                                            

use strict;

my $name;
my $seq;
my $len;
my @array;
my $min;
my $pos;

### we might need to check that the fasta sequence is all on one line.
#### need to turn this into ARGIN
open FH, "<test.fasta";
### Turn this into ARGOUT
open OUT, ">Split.fasta";
while (<FH>)
{
    if (/^>(\S+)/)
    {
        $name=$1;
    }
    if (/(\S+)/  && ! /^>/)
    {
        $seq=$1;
        $len = length($seq);
        $min = $len - 100;
        @array = split(//, $seq);
        print OUT ">$name\n$seq\n";
        print OUT ">$name\_first\n";
        for (1..99)
        {
            $pos=$_;
            print OUT "$array[$pos]";
        }
        print OUT "\n>$name\_last\n";
        for (($len-100)..$len)
        {
            $pos=$_;
            print OUT "$array[$pos]";
        }
        print OUT "\n";
    }
}
print "$len\n";
