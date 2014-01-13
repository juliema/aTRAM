use strict;
use File::Temp qw/ tempfile /;
use FindBin;
use lib "$FindBin::Bin/lib";

our $debug = 0;
our %assemblers = {};

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
	my $log_fh = shift;

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
    	my $item = pop @_[0];
        waitpid $item, 0;
    }
    return;
}

sub system_call {
	my $cmd = shift;
	my $log_fh = shift;

	unless ($log_fh) {
		$log_fh = &STDOUT;
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
	}
}

sub set_debug {
	my $debug_new = shift;
	$debug = $debug_new;
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

	close (fileIN);
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

	foreach my $contig (keys $raw_hit_matrix) {
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

sub count_partial_libraries {
	my $libname = shift;

	my $num = 0;
	while (-e "$libname.$num.1.fasta") {
		$num++;
	}
	return $num;
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

		foreach my $contig (keys $contigs) {
			my $end = $leftlength + $gaplength;
			my $start = $leftlength + 1;
			my ($startseq, $regionseq, $endseq) = split_seq ($contigs->{$contig}, $start, $end);
			contigs->{$contig} = "$startseq$endseq";
		}

		$refseq = "$left$remainder";
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
