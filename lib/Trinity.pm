#!/usr/bin/perl
use strict;
use File::Temp qw/ tempfile /;
use Module::Load;
load Assembler;

# Assembler modules need to know:
	# where to find the short reads (pass this in as a file name)
	# what the assembly parameters are. (pass this in as a hash)
# Assembler modules should return a file name for the resulting contigs.

package Trinity;

sub assembler {
	my $self = shift;
	my $short_read_file = shift;
	my $params = shift;
	my $log_fh = shift;

	my $jm = "1G";

	my ($saveout, $saveerr);
	open $saveout, ">&STDOUT";
	open $saveerr, ">&STDERR";
	open STDOUT, '>', File::Spec->devnull();
	open STDERR, '>', File::Spec->devnull();

	my $path = Assembler->find_bin("Trinity.pl");
	if ($path eq "") {
		die "couldn't find Trinity.pl ";
	}

	my ($kmer, $tempdir, $longreads, $ins_length, $exp_cov, $min_contig_len) = 0;
	if ((ref $params) =~ /HASH/) {
        if (exists $params->{"jm"}) {
			$jm = $params->{"jm"};
		}
		if (exists $params->{"tempdir"}) {
			$tempdir = $params->{"tempdir"};
		}
		if (exists $params->{"longreads"}) {
			$longreads = $params->{"longreads"};
		}
	}
	# using Trinity.pl
	print "\tassembling with Trinity...\n";
# perl ~/packages/trinityrnaseq_r20131110/Trinity.pl --seqType fa --single Pop_delt_psbA_atpA.1.blast.fasta --run_as_paired --JM 10G
	Assembler->system_call ("$path --seqType fa --single $short_read_file --run_as_paired --JM $jm --output $tempdir", $log_fh);

	open STDOUT, ">&", $saveout;
	open STDERR, ">&", $saveerr;

	return "$tempdir/Trinity.fasta";
}

sub rename_contigs {
	my $self = shift;
	my $contigfile = shift;
	my $renamefile = shift;
	my $prefix = shift;

	if ($prefix) {
		$prefix = "$prefix.";
	} else {
		$prefix = "";
	}

	open FH, "<", $contigfile;
	open OUTFH, ">", $renamefile;
	while (my $line = readline FH) {
		if ($line =~ /^>/) {
			# >comp0_c0_seq1 len=716 path=[1070:0-715]
			# >comp0_c1_seq1 len=3433 path=[4485:0-3432]
			# >comp1_c0_seq1 len=1572 path=[2936:0-1571]
			# >comp2_c0_seq1 len=4637 path=[7895:0-4636]
			$line =~ s/^>comp(\d+)_(c\d+)_seq. len=(\d+).*$/>$prefix$1$2_len_$3/;
		}
		print OUTFH $line;
	}
}

return 1;
