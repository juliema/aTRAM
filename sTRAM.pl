#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use FindBin;
use lib "$FindBin::Bin/lib";
require Subfunctions;

my $debug = 0;
if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

#required: input file names
my $short_read_archive = 0;
my $search_fasta = 0;

#optional parameters:
my $help = 0;
my $log_file = "";
my $output_file = "";
my $use_ends = 0;
my $processes = 0;
my $type = "";
my $blast_name = 0;
my $complete = 0;
my $velvet = 0;
my $protflag = 0;
my $bitscore = 70;
my $contiglength = 100;

#parameters with modifiable default values
my $ins_length = 300;
my $iterations = 5;
my $start_iter = 1;
my $exp_cov = 30;
my $evalue = 10e-10;
my $max_target_seqs = 100000000;

GetOptions ('reads=s' => \$short_read_archive,
            'target=s' => \$search_fasta,
            'start_iteration=i' => \$start_iter,
            'iterations=i' => \$iterations,
            'processes=i' => \$processes,
            'log_file=s' => \$log_file,
            'output=s' => \$output_file,
            'type=s' => \$type,
            'protein' => \$protflag,
            'blast=s' => \$blast_name,
            'use_ends' => \$use_ends,
            'debug' => \$debug,
            'velvet' => \$velvet,
            'complete' => \$complete,
            'insert_length|ins_length=i' => \$ins_length,
            'exp_coverage|expected_coverage=i' => \$exp_cov,
            'bitscore=i' => \$bitscore,
            'max_target_seqs=i' => \$max_target_seqs,
            'evalue=f' => \$evalue,
            'length=i' => \$contiglength,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

unless($short_read_archive and $search_fasta) {
    pod2usage(-msg => "Must specify a short read archive (that has already been prepared with makelibrary.pl) and a target gene in fasta form.");
}

# check to make sure that the specified short read archive exists:
if ($processes > 0) {
	unless ((-e "$short_read_archive.1.1.fasta") && (-e "$short_read_archive.1.2.fasta")) {
		pod2usage(-msg => "Short read archive does not seem to be in the format made by makelibrary.pl. Did you specify the name correctly?");
	}
	unless ((-e "$short_read_archive.$processes.1.fasta") && (-e "$short_read_archive.$processes.2.fasta")) {
		pod2usage(-msg => "Short read archives were not prepared to handle $processes processes. Try a smaller value for -processes.");
	}
} else {
	unless ((-e "$short_read_archive.1.fasta") && (-e "$short_read_archive.2.fasta")) {
		pod2usage(-msg => "Short read archive does not seem to be in the format made by makelibrary.pl. Did you specify the name correctly?");
	}
}

if ($output_file eq "") {
    $output_file = $short_read_archive;
}
if ($log_file eq "") {
	$log_file = "$output_file.log";
}

my $log_fh;
open $log_fh, ">", $log_file or die "couldn't open $log_file\n";

print $runline;
print $log_fh $runline;

my $executing_path = dirname(__FILE__);
my $cmd;


open CONTIGS_FH, ">", "$output_file.all.fasta";
truncate CONTIGS_FH, 0;
close CONTIGS_FH;

my $start_seq = "";
my $end_seq = "";
my @hit_matrices = ();
my @complete_contigs = ();

# process the target sequence file to look for the start seq and end seq.
my @target_seqs = ();

open SEARCH_FH, "<", $search_fasta;
my $name = "";
my $seq = "";
my $fullseq = "";
while (my $line=readline SEARCH_FH) {
	chomp $line;
	if ($line =~ />(.*)/) {
		if ($name ne "") {
			push @target_seqs, "$name,$seq";
		}
		$name = $1;
		$seq = "";
	} else {
		$seq .= $line;
		$fullseq .= $seq;
	}
}
push @target_seqs, "$name,$seq";

close SEARCH_FH;

my $protein = 0;
if ($protflag != 0) {
	$protein = 1;
} else {
	# if the type of target isn't specified, let's figure it out.
	if ($type ne "") {
		if ($type !~ /nucl|dna/i) {
			$protein = 1;
		}
	} else {
		$protein = is_protein($fullseq);
	}
}

my $len = 100;
if ($protein==1) {
	$len = 30;
}

# check first target seq, see if we need to make just the first bit separated for checking completeness.
my $firstseq = shift @target_seqs;
if ($firstseq =~ /(.*),(.{$len})(.*)/) {
	# split into two smaller seqs:
	$start_seq = "sTRAM_target_start";
	unshift @target_seqs, "$1,$3";
	unshift @target_seqs, "$start_seq,$2";
} else {
	# it's short enough to just go as-is
	$firstseq =~ /(.*),/;
	$start_seq = $1;
	unshift @target_seqs, $firstseq;
}

# check last target seq, see if we need to make just the first bit separated for checking completeness.
my $lastseq = pop @target_seqs;
if ($lastseq =~ /(.*),(.*)(.{$len})/) {
	# split into two smaller seqs:
	$end_seq = "sTRAM_target_end";
	push @target_seqs, "$1,$2";
	push @target_seqs, "$end_seq,$3";
} else {
	# it's short enough to just go as-is
	$lastseq =~ /(.*),/;
	$end_seq = $1;
	push @target_seqs, $lastseq;
}

# okay, let's re-assemble the file for the target fasta.
my ($TARGET_FH, $target_fasta) = tempfile(UNLINK => 1);
my @targets = ();
foreach my $line (@target_seqs) {
	$line =~ /(.*),(.*)/;
	print $TARGET_FH ">$1\n$2\n";
	push @targets, $1;
}

# make a database from the target so that we can compare contigs to the target.
my (undef, $targetdb) = tempfile(UNLINK => 1);
if ($protein == 1) {
	system_call ("makeblastdb -in $target_fasta -dbtype prot -out $targetdb.db -input_type fasta");
} else {
	system_call ("makeblastdb -in $target_fasta -dbtype nucl -out $targetdb.db -input_type fasta");
}

if ($start_iter > 1) {
	my $x = $start_iter-1;
	$search_fasta = "$output_file.".($start_iter-1).".contigs.fa";
}

# writing the header line for the results file
open RESULTS_FH, ">", "$output_file.results.txt";
print RESULTS_FH "contig\t" . join ("\t",@targets) . "\ttotal score\n";
close RESULTS_FH;

# STARTING ITERATIONS:
for (my $i=$start_iter; $i<=$iterations; $i++) {
	print ("interation $i starting...\n");
	print $log_fh ("interation $i starting...\n");

	my $hit_matrix = {};
	push @hit_matrices, \$hit_matrix;

	if ($velvet==0) {
		my $sra = "$short_read_archive";
		my $current_partial_file = "$output_file.blast.$i";
		my @partialfiles = ();
		my $p = 0;
		while ($p <= $processes) {
			$p++;
			if ($processes > 0) {
				# we are multiprocessing; make sure to name the $sra to the correct partial library
				$sra = "$short_read_archive.$p";
			}
			$current_partial_file = "$output_file.blast.$i.$p";
			push @partialfiles, $current_partial_file;

			# 1. blast to find any short reads that match the target.
			print "\tblasting short reads...\n";
			if (($protein == 1) && ($i == 1)) {
				system_call ("tblastn -max_target_seqs $max_target_seqs -db $sra.db -query $search_fasta -outfmt '6 qseqid sseqid evalue bitscore' -num_threads 8 -out $current_partial_file");
			} else {
				system_call ("blastn -task blastn -evalue $evalue -max_target_seqs $max_target_seqs -db $sra.db -query $search_fasta -outfmt '6 qseqid sseqid evalue bitscore' -num_threads 8 -out $current_partial_file");
			}

			# 2 and 3. find the paired end of all of the blast hits.
			print "\tgetting paired ends...\n";
			system_call("perl $executing_path/lib/pairedsequenceretrieval.pl $sra.#.fasta $current_partial_file $current_partial_file.fasta");
			if ($p == $processes) {
				last;
			}
		}

		# now we need to piece all of the partial files back together.
		my $fastafiles = join (" ", map {$_ . ".fasta"} @partialfiles);
		system_call("cat $fastafiles > $output_file.blast.$i.fasta");
		system_call("rm $fastafiles");

		# remove intermediate blast results
		my $readfiles = join (" ", @partialfiles);
		system_call("rm $readfiles");

		# did we not find any reads? Go ahead and quit.
		if ((-s "$output_file.blast.$i.fasta") == 0) {
			print "No similar reads were found.\n";
			last;
		}
	}

	# 4. assemble the reads using velvet
	print "\tassembling with Velvet...\n";
	system_call("velveth $output_file.velvet 31 -fasta -shortPaired $output_file.blast.$i.fasta");
	system_call("velvetg $output_file.velvet -ins_length $ins_length -exp_cov $exp_cov -min_contig_lgth 200");

	# 5. now we filter out the contigs to look for just the best ones.
	print "\tfiltering contigs...\n";

	my $blast_file = "";
	unless ($blast_name) {
		(undef, $blast_file) = tempfile(UNLINK => 1);
	} else {
		$blast_file = "$blast_name.$i";
	}

	if ($protein == 1) {
		system_call ("blastx -db $targetdb.db -query $output_file.velvet/contigs.fa -out $blast_file -outfmt '6 qseqid sseqid bitscore qstart qend sstart send qlen'");
	} else {
		system_call ("tblastx -db $targetdb.db -query $output_file.velvet/contigs.fa -out $blast_file -outfmt '6 qseqid sseqid bitscore qstart qend sstart send qlen'");
	}
	# 6. we want to keep the contigs that have a bitscore higher than $bitscore.
	open BLAST_FH, "<", $blast_file;
	while (my $line = readline BLAST_FH) {
		my ($contig, $baitseq, $score, $qstart, $qend, $sstart, $send, $qlen, undef) = split(/\s+/,$line);
		my $strand = 1;
		if ($qend > $qstart) {
			$strand = -1;
		}
		my $currscore = $hit_matrix->{$contig}->{$baitseq};
		if ($score =~ /(\d+\.\d\d)/) {
			$score = $1;
		}
		$hit_matrix->{$contig}->{"length"} = $qlen;
		if ($currscore == undef) {
			$hit_matrix->{$contig}->{$baitseq} = $strand * $score;
		} else {
			if ($currscore < $score) {
			$hit_matrix->{$contig}->{$baitseq} = $strand * $score;
			}
		}
	}
	close BLAST_FH;
	my $high_score = 0;

	my $new_matrix = {};

	foreach my $contig (keys $hit_matrix) {
		my $contig_high_score = 0;
		my $total = 0;
		$hit_matrix->{$contig}->{"strand"} = 1;
		foreach my $baitseq (@targets) {
			my $score = abs($hit_matrix->{$contig}->{$baitseq});
			$total += $score;
			if ($score > $contig_high_score) {
				$contig_high_score = $score;
				$hit_matrix->{$contig}->{"strand"} = ($hit_matrix->{$contig}->{$baitseq})/$score;
				$hit_matrix->{$contig}->{$baitseq} = $score;
			}
		}
		$hit_matrix->{$contig}->{"total"} = $total;
		if ($total > $high_score) {
			$high_score = $total;
		}
		if (($total >= $bitscore) && ($hit_matrix->{$contig}->{"length"} >= $contiglength)) {
			$new_matrix->{$contig} = $hit_matrix->{$contig};
		}
	}

	# SHUTDOWN CHECK:
	if ((keys $new_matrix) == 0) {
		print ("No contigs had a bitscore greater than $bitscore and longer than $contiglength in iteration $i: the highest bitscore this time was $high_score.\n");
		last;
	}
	$hit_matrix = $new_matrix;

	my ($SORT_FH, $sort_results) = tempfile(UNLINK => 1);
	foreach my $contig (keys $hit_matrix ) {
		print $SORT_FH $contig . "\n";
	}
	close $SORT_FH;

	# we'll use the resulting contigs as the query for the next iteration.
	$search_fasta = "$output_file.$i.contigs.fa";

	system_call("perl $executing_path/lib/findsequences.pl $output_file.velvet/contigs.fa $sort_results $search_fasta");

	# revcomping contigs with negative strand directions:
	my @contigs = ();

	open SEARCH_FH, "<", $search_fasta;
	my $name = "";
	my $seq = "";
	while (my $line=readline SEARCH_FH) {
		if ($line =~ />(.*)/) {
			if ($name ne "") {
				if ($hit_matrix->{$name}->{"strand"} < 0) {
					$seq = reverse_complement($seq);
				}
				push @contigs, ">$i"."_$name\n$seq";
			}
			$name = $1;
			chomp $name;
			$seq = "";
		} else {
			$seq .= $line;
			chomp $seq;
		}
	}
	if ($hit_matrix->{$name}->{"strand"} < 0) {
		$seq = reverse_complement($seq);
	}
	push @contigs, ">$i"."_$name\n$seq";
	close SEARCH_FH;

	# save off these resulting contigs to the ongoing contigs file.
	open CONTIGS_FH, ">>", "$output_file.all.fasta";
	print CONTIGS_FH join("\n",@contigs) . "\n";
	close CONTIGS_FH;

	open RESULTS_FH, ">>", "$output_file.results.txt";
	foreach my $contig (keys $hit_matrix) {
		my $contigname = "".($i)."_$contig";
		print RESULTS_FH "$contigname\t";
		foreach my $target (@targets) {
			if ($hit_matrix->{$contig}->{$target} == undef) {
				print RESULTS_FH "-\t";
			} else {
				my $score = abs($hit_matrix->{$contig}->{$target});
				print RESULTS_FH "$score\t";
			}
		}
		my $total = $hit_matrix->{$contig}->{"total"};
		print RESULTS_FH "$total\n";
		if ((abs($hit_matrix->{$contig}->{$start_seq}) > 70) && (abs($hit_matrix->{$contig}->{$end_seq}) > 70) && ($hit_matrix->{$contig}->{"total"} > $bitscore)) {
			push @complete_contigs, $contigname;
		}
	}
	close RESULTS_FH;

	# SHUTDOWN CHECK:
	if ($complete == 1) {
		if (@complete_contigs > 0) {
			print ("Contigs that cover the entire target sequence were found; quitting at iteration $i.\n");
			last;
		}
	}

	# if we flagged to use just the ends of the contigs, clean that up.
	if ($use_ends != 0) {
		system_call("bash $executing_path/lib/sort_contigs.sh $search_fasta");

		open FH, "<", "$search_fasta.sorted.tab";
		my @contigs = <FH>;
		close FH;
		my @new_contigs = ();

		for (my $j=0; $j<@contigs; $j++) {
			my ($len,$name,$seq) = split (/\t/,$contigs[$j]);
			if ($len > ($ins_length * 5)) {
				my $name_start = $name. "_start";
				my $name_end = $name . "_end";
				$seq =~ /^(.{$ins_length}).*(.{$ins_length})$/;
				my $seq_start = $1;
				my $seq_end = $2;
				push @new_contigs, "$ins_length\t$name_start\t$seq_start";
				push @new_contigs, "$ins_length\t$name_end\t$seq_end";
			} else {
				push @new_contigs, $contigs[$j];
			}
		}

		open FH, ">", "$search_fasta";
		foreach my $line (@new_contigs) {
			chomp $line;
			$line =~ /(.*?)\t(.*?)\t(.*)/;
			print FH ">$2\n$3\n";
		}
		close FH;
	}
}

print "\n\nContigs by target coverage:\n";

system ("cat $output_file.results.txt");
print "\nContigs containing the entire target sequence:\n\t" . join("\n\t", @complete_contigs) . "\n";
print $log_fh "Contig coverage in $output_file.results.txt\n";
print $log_fh "\nContigs containing the entire target sequence:\n\t" . join("\n\t", @complete_contigs) . "\n";

close $log_fh;

sub reverse_complement {
	my $charstr = shift;

	# reverse the DNA sequence
	my $revcomp = reverse($charstr);

	# complement the reversed DNA sequence
	$revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
	return $revcomp;
}

sub is_protein {
	my $sequence = shift;
	if ($sequence =~ /[EFILPQ]/) {
		return 1;
	} else {
		return 0;
	}
}

__END__

=head1 NAME

sTRAM.pl

=head1 SYNOPSIS

sTRAM.pl -reads shortreadfile -target target.fasta [-iterations int] [-start_iteration int] [-log_file filename] [-use_ends] [-output filename] -options

sTRAM does targeted denovo assembly of short reads to find homologs or paralogs of a target sequence.

=head1 OPTIONS
  pipeline parameters:
  -reads:     		short read archive (already run through makelibrary.pl).
  -target:          fasta file with sequences of interest.
  -output:	        optional: the prefix for the pipeline's output files (default name is the same as -reads).
  -log_file:        optional: a file to store output of the pipeline.
  -iterations:      optional: the number of pipeline iterations (default 5).
  -start_iteration: optional: if resuming from previous run, which iteration number to start from (default 0).
  -blast:			optional: saves the blast hits for the raw read library to the specified file.

  target fasta type: sTRAM, by default, tries to guess the type of fasta file specified with -target (either protein or dna).
  override options:
  -protein:         if the target sequence is a protein fasta file.
  -type:            one of the following values: dna, nucl, aa, protein.

  optional parameters:
  -use_ends:        if this flag is present, use the first and last $ins_length of long contigs in the search.
  -velvet:          if specified, skips the blast steps and performs only velvet assembly on previously-generated intermediate files.
  -complete:        if specified, automatically quits when a complete homolog is recovered.

  optional velvet assembly parameters:
  -ins_length:	    the size of the fragments used in the short-read library (default 300).
  -exp_coverage:    the expected coverage of the region for velvetg (default 30).

  optional values for searching short reads:
  -max_target_seqs: default value is 10e8.
  -evalue:          default value is 10e-10.

  optional values for blast-filtering contigs:
  -bitscore:        default value is 70.
  -length:          default value is 100.

=cut

