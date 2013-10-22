#!usr/bin/perl

$path = '/Users/juliema/Documents/Projects/TRAMiEvoBio/5-MUSCLE/*.toalign.fas';

@array =glob($path);

for $file(@array)
{
#    $count++;
#    if ($count < 5)
#    {
#  print "$file\n";
	$file =~ s/\S+MUSCLE\/(\S+).toalign.fas/$1/g;
#	print "$file\n";
	if (! exists $filehash{$file})
	{
	    %filehash=();
	    $filehash{$file}=1;
	    $numfiles++;
	    print "$numfiles\n";
	}
	system "./muscle3.8.31_i86darwin64 -in $file.toalign.fas -gapopen -2000 -out $file.muscle.out";
#    }
}
