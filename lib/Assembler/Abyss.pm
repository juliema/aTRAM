#!/usr/bin/env perl
package Abyss;
use strict;
use System;
use Parsing;
use Configuration;

my $binaries = {'abyss-pe' => "abyss-pe"};

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
#        print "spliting short read file into two files. $short_read_file $short_read_file\_1.fasta $short_read_file\_2.fasta\n";
	open FH, "<$short_read_file";
	open RD1, ">$short_read_file\_1.fasta";
	open RD2, ">$short_read_file\_2.fasta";
	my $flag=0;
	my $lib='library';
	while (<FH>) {
	    if ($flag == 1 && ! />/) { 
		my $seq=$_;
		chomp $seq;
		print RD1 "$seq\n"; 
	    }
	    if ($flag == 2 && ! />/) {
		my $seq=$_;
		chomp $seq;
		print RD2 "$seq\n"; 
	    }
	    if (/>(.*?)\/1/) {
		my $name=$1;
		$flag=1;
		print RD1;
	    }
	    if (/>(.*?)\/2/) {
		$flag=2;
		my $name=$2;
		print RD2;
	    }
	}
	Configuration::initialize();
	my ($kmer, $kmer2, $tempdir, $output_file) = 0;
	my $longreads = "";
	if ((ref $params) =~ /HASH/) {
	    if (exists $params->{"kmer"}) {
		$kmer = $params->{"kmer"};
	    }
	    if (exists $params->{"tempdir"}) {
		$tempdir = $params->{"tempdir"};
	    }
	    if (exists $params->{"longreads"}) {
		$longreads = $params->{"longreads"};
	    }
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
	    $kmer2= 2*$kmer-10;
	    ## abyss single endd
	    my $string ="v=-v k=$kmer name=$short_read_file\_temp se='$short_read_file\_1.fasta $short_read_file\_2.fasta $longreads'";
	    ### abyss paired end
##	    my $string ="v=-v k=$kmer name=$short_read_file\_temp lib='LIB' LIB='$short_read_file\_1.fasta $short_read_file\_2.fasta' se='$longreads'";
            run_command (get_bin($binaries->{'abyss-pe'}), $string,1);
	} else {
#	    my $string="v=-v k=$kmer name=$short_read_file\_temp  lib='LIB' LIB='$short_read_file\_1.fasta $short_read_file\_2.fasta'";
	    ## abyss single end
	    my $string="v=-v k=$kmer name=$short_read_file\_temp  se='$short_read_file\_1.fasta $short_read_file\_2.fasta'";
	    print "$string\n\n";
	    run_command (get_bin($binaries->{'abyss-pe'}), $string,1);
	}
	##single end abyss
	my $str = "$short_read_file\_temp-scaffolds.fa";
	my ($contigs, undef) = parsefasta ($str);
	### paired end abyss
#	my $str = "$short_read_file\_temp-contigs.fa";
#	my ($contigs, undef) = parsefasta ($str);
	print "$str\n\n";
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
#	    $contigname =~ s/^(\S+)-contigs.fa/$1/;
#	    $contigname =~ s/^>(.*)/$1/;
	    print OUTFH ">$contigname\n$sequence\n";
	}
	###### remove temp files from abyss
	`rm $short_read_file\_temp*`; 
	`rm $short_read_file\_temp*`; 
#	`rm $short_read_file_1.fasta`;
#	`rm $short_read_file_2.fasta`;
	###abyss paired end
#	`rm *.dist`;
#	`rm *.hist`;
	close OUTFH;	
	return $contigs;
}
return 1;

