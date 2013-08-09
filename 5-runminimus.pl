#!usr/bin/perl

### Need to make the output name shorter

$path = '/Applications/amos-3.1.0/bin';
$pathNames= '/Users/juliema/Documents/Projects/TRAMiEvoBio/2-RUNAMOS/*.FASTAOUT.fas';

@array = glob("$pathNames");

for $file(@array)
{
    $count++;
    print "$count\n";
    system "$path/toAmos -s $file -o $file.afg";
    system "$path/minimus $file.afg"; 
}

print "There are $count files\n";

