#!/usr/bin/env perl
package Subfunctions;
use strict;
use File::Temp qw/ tempfile /;
use FindBin;
use lib "$FindBin::Bin/lib";

BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw( $log_fh parse_config find_bin timestamp exit_with_msg fork_cmd wait_for_forks printlog system_call debug set_debug set_log parsefasta sortfasta flattenfasta make_hit_matrix process_hit_matrix set_multiplier get_multiplier map_to_shard set_total_shards get_total_shards get_max_shard is_protein percentcoverage split_seq );
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw($total_shards %assemblers $multiplier);
}

our $debug = 0;
our %assemblers = {};
our $log_fh = 0;
our $total_shards = 0;
our @primes = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151);
our $multiplier = 0;


sub parse_config {
	open FH, "<", "$FindBin::Bin/config.txt";
	foreach my $line (<FH>) {
		$line =~ s/(#.*)$//;
		if ($line =~ /(.*)=(.*)$/) {
			my $name = $1;
			my $path = $2;
			$assemblers{$name} = "$path";
		}
	}
	close FH;
}

sub find_bin {
	my $cmd = shift;

	if (exists $assemblers{$cmd}) {
		print "found $cmd: at $assemblers{$cmd}\n";
	}
}

sub timestamp {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $mon++;
    $mon = sprintf("%02d", $mon);
    $min = sprintf("%02d", $min);
    $sec = sprintf("%02d", $sec);
    $hour = sprintf("%02d", $hour);
    $mday = sprintf("%02d", $mday);

    $year -= 100;
    my $time = "$hour:$min:$sec";
    my $date = "$year$mon$mday";
    return "$date $time";
}

sub exit_with_msg {
	my $msg = shift;
	print STDERR "$msg\n";
	exit 1;
}

sub fork_cmd {
	my $cmd = shift;
	print $log_fh ("\t$cmd\n");
    my $child_pid = fork();
    unless ($child_pid) { #child process
		exec ($cmd);
    } else { #parent process
        return $child_pid;
    }
}

sub wait_for_forks {
    while (@{@_[0]} > 0) {
     	debug ("waiting for " . join(", ", @{@_[0]}) . "\n");
		my $item = pop @{@_[0]};
        waitpid $item, 0;
    }
    return;
}

sub printlog {
	my $msg = shift;

	$msg = timestamp() . ": " . $msg . "\n";
	print $msg;
	if ($log_fh) {
        select($log_fh);
        $|++;
		print $log_fh $msg;
		select(STDOUT);
	}
}

sub system_call {
	my $cmd = shift;

	if ($log_fh == 0) {
		open my $std_log, ">&", STDOUT;
		$log_fh = $std_log;
	}

	print $log_fh ("\t$cmd\n");
	my ($saveout, $saveerr);
	if ($debug == 0) {
		open $saveout, ">&STDOUT";
		open $saveerr, ">&STDERR";
		open STDOUT, '>', File::Spec->devnull();
		open STDERR, '>', File::Spec->devnull();
	}
	my $exit_val = eval {
		system ($cmd);
	};

	if ($debug == 0) {
		open STDOUT, ">&", $saveout;
		open STDERR, ">&", $saveerr;
	}

	if ($exit_val != 0) {
		print "System call \"$cmd\" exited with $exit_val\n";
		exit;
	}
	return $exit_val;
}

sub debug {
	my $msg = shift;
	if ($debug) {
		print "$msg";
		print $log_fh ("DEBUG: $msg");
	}
}

sub set_debug {
	my $debug_new = shift;
	$debug = $debug_new;
}

sub set_log {
	my $log_new = shift;
	$log_fh = $log_new;
}

sub parsefasta {
	my $fastafile = shift;

	my $taxa = {};
	my @taxanames = ();
	open fileIN, "<", "$fastafile";
	my $input = readline fileIN;
	my $taxonlabel = "";
	my $sequence = "";
	while ($input !~ /^\s*$/) {
		if ($input =~ /^>(.+)\s*$/) {
			$taxonlabel = $1;
			$taxonlabel =~ s/\s+/_/g;
			push @taxanames, $taxonlabel;
		} else {
			$input =~ /^\s*(.+)\s*$/;
			$taxa->{$taxonlabel} .= $1;
		}
		$input = readline fileIN;
	}

	close fileIN;
	return $taxa, \@taxanames;
}

sub sortfasta {
	my $fastafile = shift;
	my $outfile = shift;
	my $separator = shift;

	unless ($separator) {
		$separator = '\n';
	}

	my ($taxa, $taxanames) = parsefasta ($fastafile);
	my @sortedtaxa = sort @$taxanames;
	open fileOUT, ">", "$outfile";
	foreach my $taxon (@sortedtaxa) {
		print fileOUT ">$taxon$separator$taxa->{$taxon}\n";
	}
	close fileOUT;
}

sub flattenfasta {
	my $fastafile = shift;
	my $outfile = shift;
	my $separator = shift;

	unless ($separator) {
		$separator = '\n';
	}

	my ($taxa, $taxanames) = parsefasta ($fastafile);
	open fileOUT, ">", "$outfile";
	foreach my $taxon (@$taxanames) {
		print fileOUT ">$taxon$separator$taxa->{$taxon}\n";
	}
	close fileOUT;
}

sub make_hit_matrix {
	my $blast_file = shift;
	my $hit_matrix = shift;
	my $iter = shift;

	open BLAST_FH, "<", $blast_file;
	while (my $line = readline BLAST_FH) {
		my ($contig, $baitseq, $score, $qstart, $qend, $sstart, $send, $qlen, undef) = split(/\s+/,$line);
		my $strand = 1;
		if ((($qend-$qstart) / ($send-$sstart)) < 0) {
			$strand = -1;
		}
		my $currscore = $hit_matrix->{$contig}->{$baitseq};
		if ($score =~ /(\d+\.\d\d)/) {
			$score = $1;
		}
		$hit_matrix->{$contig}->{"length"} = $qlen;
		$hit_matrix->{$contig}->{"iteration"} = $iter;
		if ($currscore == undef) {
			$hit_matrix->{$contig}->{$baitseq} = $strand * $score;
		} else {
			if (abs($currscore) < $score) {
			$hit_matrix->{$contig}->{$baitseq} = $strand * $score;
			}
		}
	}
	close BLAST_FH;
}

sub process_hit_matrix {
	my $raw_hit_matrix = shift;
	my $targets = shift;
	my $bitscore = shift;
	my $contiglength = shift;
	my $hit_matrix = shift;
	my $contig_names = shift;

	my $high_score = 0;

	# clean up the hit matrix: only keep hits that meet the bitscore threshold.

	foreach my $contig (keys %$raw_hit_matrix) {
		my $contig_high_score = 0;
		my $total = 0;
		$raw_hit_matrix->{$contig}->{"strand"} = 1;
		foreach my $baitseq (@$targets) {
			my $partscore = abs($raw_hit_matrix->{$contig}->{$baitseq});
			if ($partscore) {
				my $partstrand = ($raw_hit_matrix->{$contig}->{$baitseq})/$partscore;
				if ($partscore > 0) {
					# separate out the score and the strand for this part:
					$raw_hit_matrix->{$contig}->{$baitseq} = $partscore;
					$total += $partscore;
					if ($partscore > $contig_high_score) {
						# if this is the best score for the contig, set the contig_high_score and set the strand to this strand.
						$contig_high_score = $partscore;
						$raw_hit_matrix->{$contig}->{"strand"} = $partstrand;
					}
				}
			}
		}
		$raw_hit_matrix->{$contig}->{"total"} = $total;
		if ($total > $high_score) {
			$high_score = $total;
		}
		if (($total >= $bitscore) && ($raw_hit_matrix->{$contig}->{"length"} >= $contiglength)) {
			# check to see if this contig is already in the hit_matrix:
			if (exists $hit_matrix->{$contig}) {
				if (($hit_matrix->{$contig}->{"length"} == $raw_hit_matrix->{$contig}->{"length"}) && ($hit_matrix->{$contig}->{"total"} == $raw_hit_matrix->{$contig}->{"total"})) {
					# it's the same, don't add it.
					next;
				}
			}
			$hit_matrix->{$contig} = $raw_hit_matrix->{$contig};
			push @$contig_names, $contig;
		}
	}
	return $high_score;
}

#### MAPREDUCE FUNCTIONS

sub set_multiplier {
	my $factor = shift;
	$multiplier = $primes[$factor % (@primes)];
	return $multiplier;
}

sub get_multiplier {
	if ($multiplier != 0) {
		$multiplier = $primes[0];
	}
	return $multiplier;
}

sub map_to_shard {
	my $name = shift;

	$name =~ s/\/\d//;
	$name =~ s/#.+$//;
	$name =~ tr/0-9//dc;

	$name =~ /.*(\d{8})$/;
	$name = $1 * get_multiplier();

	return $name % $total_shards;
}

sub set_total_shards {
	$total_shards = shift;
	return $total_shards;
}

sub get_total_shards {
	my $dbname = shift;

	if ($total_shards == 0) {
		my $num = 0;
		while (-e "$dbname.$num.1.fasta") {
			$num++;
		}
		$total_shards = $num;
	}
	return $total_shards;
}

sub get_max_shard {
	return get_total_shards() - 1;
}

sub is_protein {
	my $sequence = shift;
	if ($sequence =~ /[EFILPQ]/) {
		return 1;
	} else {
		return 0;
	}
}

sub percentcoverage {
	my $reffile = shift;
	my $contigfile = shift;
	my $outname = shift;
	my $aligner = shift;

	my (undef, $catfile) = tempfile(UNLINK => 1);
	system ("cat $reffile $contigfile > $catfile");

	if ($aligner eq "mafft") {
		system ("mafft  $catfile > $outname.align.fasta");
	} else {
		system ("muscle -in $catfile -out $outname.align.fasta");
	}

	open REF_FH, "<", $reffile;
	my $ref = readline REF_FH;
	$ref =~ />(.+)$/;
	my $refname = $1;
	close REF_FH;

	# parse the output file: save the reference as a separate sequence, put the others into an array.
	my ($contigs, $taxa) = parsefasta ("$outname.align.fasta");
	my $refseq = delete $contigs->{"$refname"};

	# as long as there are still gaps in the reference sequence, keep removing the corresponding positions from the contigs.
	while ($refseq =~ /(\w*)(-+)(.*)/) {
		my $left = $1;
		my $gap = $2;
		my $remainder = $3;
		my $leftlength = length $left;
		my $gaplength = length $gap;

		foreach my $contig (keys %$contigs) {
			my $end = $leftlength + $gaplength;
			my $start = $leftlength + 1;
			my ($startseq, $regionseq, $endseq) = split_seq ($contigs->{$contig}, $start, $end);
			$contigs->{$contig} = "$startseq$endseq";
		}

		$refseq = "$left$remainder";
	}
	# align the ends of the contig seqs
	foreach my $contig (keys %$contigs) {
		if ((length $contigs->{$contig}) > (length $refseq)) {
			# if the contig seq is longer than the refseq, truncate it
			$contigs->{$contig} = substr($contigs ->{$contig}, 0, length $refseq);
		} else {
			# if the contig seq is shorter than the refseq, pad it with gaps.
			$contigs->{$contig} .= "-" x ((length $contigs->{$contig}) - (length $refseq));
		}
	}

	$contigs->{"reference"} = $refseq;
	return $contigs;
}


sub split_seq {
    my $seq = shift;
    my $start = shift;
    my $end = shift;
    my $max = 30000;
        my $seqlen = length ($seq);
        my $startseq = "";
        my $regionseq = "";
        my $endseq = "";

        my $currstart = $start-1;
        my $currend = $end;
        while ($currstart > $max) {
                $seq =~ /^(.{$max})(.*)$/;
                $startseq .= $1;
                $seq = $2;
                $currstart -= $max;
                $currend -= $max;
        }
        if ($currstart > 0) {
                $seq =~ /^(.{$currstart})(.*)$/;
                $startseq .= $1;
                $seq = $2;
                $currstart = 1;
                $currend = $end - (length ($startseq));
        }

        my $regionsize = $end - $start + 1;
        while ($regionsize > $max) {
                $seq =~ /^(.{$max})(.*)$/;
                $regionseq .= $1;
                $seq = $2;
                $currstart -= $max;
                $currend -= $max;
                $regionsize -= $max;
        }
        if ($regionsize > 0) {
                $seq =~ /^(.{$regionsize})(.*)$/;
                $regionseq .= $1;
                $endseq = $2;
        }
        return ($startseq, $regionseq, $endseq);
}


return 1;
