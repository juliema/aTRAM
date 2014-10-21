#!/usr/bin/env perl
package Abyss;
use strict;
use System;
use Parsing;
use Configuration;

# Assembler modules need to know:
	# where to find the short reads (pass this in as a file name)
	# what the assembly parameters are. (pass this in as a hash)
# Assembler modules should return a hash of the resulting contigs.

# Hash of assembler's required binaries
#my $binaries = {velveth => "velveth", velvetg => "velvetg"};
my $binaries = {abyss-pe => "abyss-pe"};


sub get_binaries {
	return $binaries;
}

sub assembler {
	my $self = shift;
	my $short_read_file = shift;
	my $params = shift;

###############
#### need to split the reads into two files return  $shortreadfile.read1  $shortreadfile.read2
##############

	Configuration::initialize();
#	my ($kmer, $tempdir, $ins_length, $exp_cov, $min_contig_len, $output_file) = 0;
	my ($k, $tempdir, $output_file) = 0;
	my $longreads = "";

	if ((ref $params) =~ /HASH/) {
	    if (exists $params->{"k"}) {
		$k= $params->{"k"};
	    }
	    else {$k = 64}
	    if (exists $params->{"tempdir"}) {
		$tempdir = $params->{"tempdir"};
	    }
	    if (exists $params->{"longreads"}) {
		$longreads = $params->{"longreads"};
	    }
#		if (exists $params->{"ins_length"}) {
#			$ins_length = $params->{"ins_length"};
#		}
#		if (exists $params->{"exp_cov"}) {
#			$exp_cov = $params->{"exp_cov"};
#		}
#		if (exists $params->{"min_contig_len"}) {
#			$min_contig_len = $params->{"min_contig_len"};
#		}
	    if (exists $params->{"output"}) {
		$output_file = $params->{"output"};
	    }
	    if (exists $params->{"log_file"}) {
		set_log($params->{"log_file"});
	    }
	}
	# using ABySS
	# truncate ABySS log file if it already exists
	truncate "$tempdir/Log", 0;
	if ($longreads ne "") {
	    #abyss-pe k=64 name=aTRAM-test lib='pe160' pe160='/path/reads.1.fasta /path/reads.2.fasta'                                               
	    run_command (get_bin($binaries->{abyss-pe}), "$k $name $lib $lib='$shortreadfile.read1  $shortreadfile.read2",1);
	} else {
	    #abyss-pe k=64 name=aTRAM-test lib='pe160' long=long1 pe160='/path/reads.1.fasta /path/reads.2.fasta' long1='/path/aTRAM-test-contigs.fa'
	    run_command (get_bin($binaries->{abyss-pe}), " $k, $name, $lib, long=$longreads $lib=$shortreadfile.read1  $shortreadfile.read2",1);
	}
	my ($contigs, undef) = parsefasta ("$tempdir/contigs.fa");
	# copy ABySS log output to logfile.
	open LOGFH, "<:crlf", "$tempdir/Log";
	printlog ("Abyss log:");
	foreach my $line (<LOGFH>) {
	    chomp $line;
	    printlog ($line);
	}
	printlog ("end Abyss log");
	close LOGFH;
	
	open OUTFH, ">", $output_file;
	foreach my $contigname (keys %$contigs) {
	    my $sequence = $contigs->{$contigname};
	    #aTRAM-test-contigs.fa
#	    $contigname =~ s/^NODE_(\d+)_length_(\d+)_cov_(\d+\.\d).*$/$1_len_$2_cov_$3/;
	    $contigname =~ s/^(\S+)-contigs.fa/$1/;
	    print OUTFH ">$contigname\n$sequence\n";
	}
	close OUTFH;
	return $contigs;
}
return 1;



#	if ($longreads ne "") {
#		run_command (get_bin($binaries->{velveth}), "$tempdir $kmer -fasta -shortPaired $short_read_file -long $longreads", 1);
#	} else {
#		run_command (get_bin($binaries->{velveth}), "$tempdir $kmer -fasta -shortPaired $short_read_file", 1);
#	}
#	run_command (get_bin($binaries->{velvetg}), "$tempdir -ins_length $ins_length -exp_cov $exp_cov -min_contig_lgth $min_contig_len", 1);
#	my ($contigs, undef) = parsefasta ("$tempdir/contigs.fa");
#
#	# copy Velvet log output to logfile.
#	open LOGFH, "<:crlf", "$tempdir/Log";
#	printlog ("Velvet log:");
#	foreach my $line (<LOGFH>) {
#		chomp $line;
#		printlog ($line);
#	}
#	printlog ("end Velvet log");
#	close LOGFH;
#
#	open OUTFH, ">", $output_file;
#	foreach my $contigname (keys %$contigs) {
#		my $sequence = $contigs->{$contigname};
#		#aTRAM-test-contigs.fa
#		#NODE_41_length_2668_cov_4.901050
#		$contigname =~ s/^NODE_(\d+)_length_(\d+)_cov_(\d+\.\d).*$/$1_len_$2_cov_$3/;
#		print OUTFH ">$contigname\n$sequence\n";
#	}
#	close OUTFH;
#	return $contigs;
#}
#return 1;


#To do the basic call with just paired end reads:
#abyss-pe k=64 name=aTRAM-test lib='pe160' pe160='/path/reads.1.fasta /path/reads.2.fasta'

#to add the contigs back:
#abyss-pe k=64 name=aTRAM-test lib='pe160' long=long1 pe160='/path/reads.1.fasta /path/reads.2.fasta' long1='/path/aTRAM-test-contigs.fa'


#explanation
#abyss-pe  ->  program call
#k=64  ->  this is to define k-mer, here I am using 64
#name=aTRAM-test  ->  this defines the file handle name, I called it aTRAM-test in the above examples
#lib='pe160' ->  this lets make a name to define a library, I called it pe160 in the above examples 
#pe160='/path/reads1 /path/reads2'  ->  this applies the data to a named library
#long='long1'  ->  this allows you to tell abyss you have long reads for scaffolding, I called it long1 in the above examples
#long1='aTRAM-test-contigs.fa'  ->  this allows you to defined your contigs from a previous run to the name "example_long"  I called it long1 in the example above

#Here is an example report from ABySS, this way you can see all the files output
#n       n:500   n:N50   min     N80     N50     N20     max     sum
#314410   105     35      503     605     795     1414    2608    87279   aTRAM-test-unitigs.fa
#13975   228     61      500     589     887     2348    5581    211832  aTRAM-test-contigs.fa
#13959   222     55      500     592     983     3510    5581    215055  aTRAM-test-scaffolds.fa

#The contigs file, will be a multi-fasta file.
