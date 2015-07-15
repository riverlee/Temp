#!/usr/bin/env perl

=head1 NAME

xhyb_check.pl

=head1 SYNOPSIS

xhyb_check.pl
  --query_file <file>
  --database_file <file>
 [--ignore_file <file>]
 [--common_file <file>]
 [--hash_length <number>]
 [--merge 0]
 [--sequence_only]
 [--db_start <number>]
 [--db_stop <number>]
 [--verbose [verbosity_level]]

=head1 DESCRIPTION

This script looks for Nmers in the query file in the database file and reports
all hits.  When a hit is found, the hit acts like a seed and bases upstream and
downstream are checked for matches and the hit is extended if possible.  You can
also specify an ignore file with a list of common sequences to ignore (ie, repeats,
LINES, etc).

=head1 FOR GEEKS

First, all the Nmers in the sequences to ignore is stored in a hash.  Next, all
the Nmers in the query sequences are stored in a hash minus the ignored Nmers.
Last, the script walks through database sequences one position at a time and
checks the query hash.  When a hit is found, the bases upstream and downstream
are compared to the query sequence and the hit is extended when possible.  Any
expected subsequent Nmer hits to the same area are ignored.

=head1 OPTIONS

=over 2

=item B<--query_file>  file

Required.  This is a file of sequences the program will store in a hash.  It is
recommended that this is the smaller of the two sets of files.  You may specify
more than one query file option.

=item B<--database_file>  file

Required.  This is a FASTA file of sequences the program will use to walk through
and check for hits in the query.  It is recommended that this is the larger of the
two sets of files.  You may specify more than one database file option.

=item B<--ignore_file>  file

Optional.  This is a file of sequences the program will use to ignore.  You
may specify more than one ignore file option.

=item B<--common_file> file

Optional.  This is a file of subsequences the program will use to ignore.  The
intention here is that common subseqs make the probes perform poorly.  The length
of the sequences in this file should be the same length as your hash.

=item B<--hash_length>  number

Optional.  Default is 25.  This is the length of the keys in the hash.  Generally, a
shorter hash length will result in more memory usage because there are more Nmers
to store.

=item B<--merge>  number

Optional.  Default is 1.  Use 0 to not merge consecutive hits and maintain individual
hits.

=item B<--allow_mismatches>  [number]

Optional.  Allows [number] mismatches when extending the hit.  If no number is given,
unlimited mismatches are allowed in the hit.

=item B<--sequence_only>

Optional.  Use this option to only output the originating probe and the alignment
subsequence.  One line will be output for each hit.  This option is normally used
to save disk space when using a hash length somewhat smaller than the length of the
probe (eg, 16mers with 25mer probes).

=item B<--pad_alignment>

Optional.  Unaligned bases upstream and downstream in the query sequence will be
represented as ".".

=item B<--use_xhyb_rules>

Optional.  Only output hits that follow Chip Design's crosshyb rules which are
23 out of 25 matches, 16-mer hit, or 12mer + 8mer hit.

=item B<--strand>  +|-

Optional.  Only look on specified strand of the target for matches.

=item B<--db_start>  number

Optional.  In the event you need to partition a large sequence into smaller parts,
you may use this option to specify a starting position in the database.  Note that
this option may be useless without specifying the --db_stop.

=item B<--db_stop>  number

Optional.  In the event you need to partition a large sequence into smaller parts,
you may use this option to specify a stopping position in the database.  Note that
this option may be useless without specifying the --db_start.

=item B<--use_qnames>

Optional.  Output the query as qnames instead of the query sequence.

=item B<--verbose>  [verbosity_level]

Optional.  If unspecified, verbosity level is defaulted to 1.  Sets the verbosity
level.

=back

=head1 COMPLAINTS

Brant Wong, brant_wong@affymetrix.com

=head1 COPYRIGHT

Affymetrix, Inc., 2006.  All rights reserved.

=cut


use strict;
use Getopt::Long;
use Neo::Fasta;
use Affx::Util::Parameter;
use Affx::SS::Tools;

my ($hash_length) = (25);
my (@query_files, @database_files, @common_files, @ignore_files);
my ($merge) = (1);
my ($allow_mismatches);
my ($sequence_only) = (0);
my ($pad_alignment, $use_xhyb_rules) = (0, 0);
my ($strand) = ('');
my ($db_start, $db_stop) = (1, 0);
my ($use_qnames) = (0);
my ($verbose);

unless (GetOptions("hash_length=i"       => \$hash_length,
		   "query_file=s"        => \@query_files,
		   "database_file=s"     => \@database_files,
		   "common_file=s"       => \@common_files,
		   "ignore_file=s"       => \@ignore_files,
		   "merge=i"             => \$merge,
		   "allow_mismatches:i"  => \$allow_mismatches,
		   "sequence_only"       => \$sequence_only,
		   "pad_alignment"       => \$pad_alignment,
		   "use_xhyb_rules"      => \$use_xhyb_rules,
		   "strand=s"            => \$strand,
		   "db_start=i"          => \$db_start,
		   "db_stop=i"           => \$db_stop,
		   "use_qnames"          => \$use_qnames,
		   "verbose:i"           => \$verbose,
		  )) {die "GetOptions error\n";}


unless ($hash_length) {die "Must supply --hash_length\n";}
unless (@query_files) {die "Must supply --query_file (smaller set)\n";}
unless (@database_files) {die "Must supply --database_file (larger set)\n";}

if (defined $verbose && not $verbose) {$verbose = 1;}
if (defined $allow_mismatches && not $allow_mismatches) {$allow_mismatches = 999999999;}


my %query_hash;
my @query_seqs;
my %seq_hash;
my %common_hash;
my %ignore_hash;
my %qnames;

warn "Starting:  " . `date`;
if (@common_files) {
  if ($verbose) {warn "Storing hash of common sequences ## " . `date`;}
  foreach my $common_file (@common_files) {
    if ($verbose > 1) {warn "Looking through $common_file ## " . `date`;}
    my $common_file_F = Affx::SS::Tools::NewFileHandle("$common_file");
    while (my $common_line = <$common_file_F>) {
      if ($common_line =~ /^>/) {next;}
      chomp $common_line;
      $common_line = uc($common_line);
      $common_hash{$common_line} = 1;
    }
  }
}


if ($verbose) {warn "Loading hash from query files ## " . `date`;}
foreach my $query_file (@query_files) {
  if ($verbose > 1) {warn "Looking through $query_file ## " . `date`;}
  my $qname = '';
  my $query_seq = '';
  my $query_file_F = Affx::SS::Tools::NewFileHandle("$query_file");
  while (my $query_line = <$query_file_F>) {
    if ($query_line =~ /^>(\S+)/) {
      if ($query_seq) {
	$query_hash{$query_seq} = 1;
	if ($use_qnames) {$qnames{$query_seq}{$qname} = 1;}
	$query_seq = '';
      }

      # Store qname to flag that this is a fasta file.
      $qname = $1;
    }
    else {
      chomp $query_line;
      $query_line = uc($query_line);
      $query_line =~ s/\s+//g;

      # If this is a fasta file, we want to combine the lines.
      if ($qname) {$query_seq .= $query_line;}

      # If this is not a fasta file, assume each line is a query.
      else {
	$query_hash{$query_line} = 1;
	if ($use_qnames) {$qnames{$query_line}{$query_line} = 1;}
      }
    }
  }
  if ($query_seq) {
    $query_hash{$query_seq} = 1;
    if ($use_qnames) {$qnames{$query_seq}{$qname} = 1;}
  }
}


if (@ignore_files && $verbose) {warn "Checking sequences in query files against ignore files ## " . `date`;}
foreach my $ignore_file (@ignore_files) {
  if ($verbose > 1) {warn "Looking through $ignore_file ## " . `date`;}
  my $ignore_file_F = Affx::SS::Tools::NewFileHandle("$ignore_file");
  while (my $ignore_line = <$ignore_file_F>) {
    if ($ignore_line =~ /^>/) {next;}
    chomp $ignore_line;
    $ignore_line = uc($ignore_line);
    if (exists $query_hash{$ignore_line}) {delete $query_hash{$ignore_line};}
  }
}

if ($verbose) {warn "Hashing subseqs in queries ## " . `date`;}
my $max_query_size = 0;
foreach my $query_seq (keys %query_hash) {
  if ($max_query_size < length($query_seq)) {$max_query_size = length($query_seq);}
  ProcessHash(seq => \$query_seq, operation => 'load');
  delete $query_hash{$query_seq};
}


if ($verbose) {warn "Checking hash with database files ## " . `date`;}
foreach my $database_file (@database_files) {
  my $database_file_F = Affx::SS::Tools::NewFileHandle("$database_file");
  my $database_file_FA = new Neo::Fasta();
  $database_file_FA->Normalize(0);
  $database_file_FA->Read($database_file_F, \&CheckHash);
}
warn "Finished:  " . `date`;


sub CheckHash {

  my ($label, $seq) = @_;

  if ($strand eq '+') {ProcessHash(label => $label, seq => $seq, operation => 'check');}
  elsif ($strand eq '-') {ProcessHash(label => $label, seq => $seq, operation => 'check', revcomp => 1);}
  else {
    ProcessHash(label => $label, seq => $seq, operation => 'check');
    ProcessHash(label => $label, seq => $seq, operation => 'check', revcomp => 1);
  }
}


sub ProcessHash {

  my ($label, $seq, $operation, $revcomp) = Affx::Util::Parameter::getparams(["label", "seq", "operation", "revcomp"], @_);

  if ($revcomp) {
    $revcomp = '-';
    $seq = Neo::Fasta::ReverseComp($seq);
  }
  else {$revcomp = '+';}

  if ($operation eq 'check' && $verbose > 1) {warn "Looking through $label on $revcomp strand ## " . `date`;}

  # Get the total length of the sequence.
  my $seq_length = length($$seq);

  # This is used to keep the location in the sequence.
  my $location = ($operation eq 'check') ? $db_start - 1 : 0;

  # Initiate hash of previous hits.
  my %previous_hits;

  # Continue as long as we have enough bases in the sequence to extract.
  while (++$location <= $seq_length - $hash_length + 1) {

    # Extract the sub sequence.
    my $subseq = substr($$seq, $location - 1, $hash_length);

    # If we are just loading the sequence hash...
    if ($operation eq 'load') {
      unless (exists $common_hash{$subseq}) {$seq_hash{$subseq}{$$seq}{$location} = 1;}
    }

    # If we are checking against the sequence hash...
    elsif ($operation eq 'check') {

      # Check to see if the subseq exists in the sequence hash.
      if (not(exists $common_hash{$subseq}) && exists $seq_hash{$subseq}) {

	# Set flag for beginning sequence.
	my $beginning = ($location + $hash_length - $max_query_size - 1 < 0) ? 1 : 0;

	# Extract context sequence to compare.
	my $context_seq_start = $beginning ? 0 : $location + $hash_length - $max_query_size - 1;
	my $context_seq = substr($$seq, $context_seq_start, $max_query_size * 2 - $hash_length);
	my @split_context_seq = split(//, $context_seq);

	# Go through the labels and locations to record a new hit.
	foreach my $probe (keys %{$seq_hash{$subseq}}) {

	  # Split probe for base matching.
	  my @split_probe = split(//, $probe);

	  foreach my $probe_location (keys %{$seq_hash{$subseq}{$probe}}) {

	    unless ($merge) {
	      if ($use_qnames) {
		unless (exists $qnames{$probe}) {die "Could not find qname for probe:  $probe\n";}
		foreach my $qname (keys %{$qnames{$probe}}) {
		  if ($sequence_only) {print join("\t", $qname, $subseq), "\n";}
		  else {print join("\t", $label, $revcomp, $location, $hash_length, $qname, $probe_location, $subseq), "\n";}
		}
	      }
	      else {
		if ($sequence_only) {print join("\t", $probe, $subseq), "\n";}
		else {print join("\t", $label, $revcomp, $location, $hash_length, $probe, $probe_location, $subseq), "\n";}
	      }
	    }

	    else {
	      if (exists $previous_hits{$probe}{$probe_location}) {next;}
	      $previous_hits{$probe}{$probe_location} = 1;

	      my $probe_match_start = $probe_location;
	      my $probe_match_stop = $probe_location + $hash_length - 1;
	      my $align_seq = $subseq;
	      my $offset = 1;

	      # First thing to do is to look as far back as possible for matches.
	      my $test_probe_location = $probe_location - 1;
	      for (my $previous_base = $beginning ? $location - 1 : $max_query_size - $hash_length; $test_probe_location > 0 && $previous_base > 0; $test_probe_location--, $previous_base--) {

		# We found a match.
		if ($split_probe[$test_probe_location - 1] eq $split_context_seq[$previous_base - 1]) {
		  $align_seq = "$split_context_seq[$previous_base - 1]$align_seq";
		  $probe_match_start -= $offset;
		  $offset = 1;
		}
		else {
		  $align_seq = ".$align_seq";
		  $offset++;
		}
	      }

	      # Pre-pad alignment for hits close to the end of the db.
	      if ($test_probe_location) {$align_seq = ('.' x $test_probe_location) . $align_seq;}

	      # Reset offset.
	      $offset = 1;

	      # Next go forward and look for matches.
	      $test_probe_location = $probe_location + $hash_length;
	      my $next_base = $beginning ? $location + $hash_length - 1 : $max_query_size;
	      for (; $test_probe_location <= length($probe) && $next_base < @split_context_seq; $test_probe_location++, $next_base++) {

		# We found a match.
		if ($split_probe[$test_probe_location - 1] eq $split_context_seq[$next_base]) {
		  $align_seq .= "$split_context_seq[$next_base]";
		  $probe_match_stop += $offset;
		  $offset = 1;
		}
		else {
		  $align_seq .= '.';
		  $offset++;
		}
	      }

	      # Post-pad alignment for hits close to the end of the db.
	      if ($test_probe_location < length($probe)) {$align_seq .= '.' x (length($probe) - $test_probe_location + 1);}
	      elsif ($test_probe_location == length($probe) && $next_base == @split_context_seq) {$align_seq .= '.';}

	      # Removed padded alignment if not requested.
	      unless ($pad_alignment) {$align_seq =~ s/^\.|\.$//g;}

	      # So we don't modify our location when revcomping.
	      my $db_location = $location - ($probe_location - $probe_match_start);

	      if ($revcomp eq '-') {

		# Explanation of formula:
		#  +1 to probe stop - probe start to get length of the alignment
		#  -1 from db_location to get the actual base position of the right side of the alignment
		#  +1 to seq_length because we want to be in 1-based index when moving to the opposite strand
		$db_location = $seq_length - ($db_location + ($probe_match_stop - $probe_match_start + 1) - 1) + 1;
	      }

	      # Skip if too many mismatches.
	      my $okay_to_print = 1;
	      if ($allow_mismatches) {
		my $num_mismatches = $align_seq =~ tr/\./\./;
		if ($num_mismatches > $allow_mismatches) {$okay_to_print = 0;}
	      }

	      # Check for xhyb rules.
	      if ($use_xhyb_rules && $okay_to_print) {

		# 23 out of 25.
		my $match_count = $align_seq =~ tr/ACGT/ACGT/;
		unless ($match_count >= 23 ||

			# 16mer
			$align_seq =~ /[^\.]{16}/ ||

			# 8 + 12
			$align_seq =~ /[^\.]{8,}.*[^\.]{12,}/ || $align_seq =~ /[^\.]{12,}.*[^\.]{8,}/) {

		  # Skip if we don't have one of these rules.
		  $okay_to_print = 0;
		}
	      }


	      if ($okay_to_print) {
		if ($use_qnames) {
		  unless (exists $qnames{$probe}) {die "Could not find qname for probe:  $probe\n";}
		  foreach my $qname (keys %{$qnames{$probe}}) {
		    if ($sequence_only) {print join("\t", $qname, $align_seq), "\n";}
		    else {print join("\t", $label, $revcomp, $db_location, $db_location + $probe_match_stop - $probe_match_start, $probe_match_stop - $probe_match_start + 1, $qname, $probe, $probe_match_start, $probe_match_start + $probe_match_stop - $probe_match_start, $align_seq), "\n";}
		  }
		}
		else {
		  if ($sequence_only) {print join("\t", $probe, $align_seq), "\n";}
		  else {print join("\t", $label, $revcomp, $db_location, $probe_match_stop - $probe_match_start + 1, $probe, $probe_match_start, $align_seq), "\n";}
		}
	      }
	    }
	  }
	}
      }
    }

    # Advance hits we've seen.
    my %new_previous_hits;
    foreach my $probe (keys %previous_hits) {
      foreach my $probe_location (keys %{$previous_hits{$probe}}) {
	if ($probe_location <= length($probe) - $hash_length) {$new_previous_hits{$probe}{$probe_location + 1} = 1;}
      }
    }
    %previous_hits = %new_previous_hits;

    # Exit loop if db_stop.
    if ($db_stop && $location >= $db_stop - $hash_length + 1) {last;}
  }
}

