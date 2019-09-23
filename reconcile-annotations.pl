#! /usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Carp::Always;
use Data::Dumper;

use Set::IntervalTree;

# use FindBin;
# use lib "$FindBin::Bin";
# use Xyzzy;

use constant { TRUE => 1, FALSE => 0 };

sub max {
  my ($x,$y) = @_;
  if ( !defined($x) ) {
    return $y;
  } elsif ( !defined($y) ) {
    return $x;
  } elsif ( $x > $y ) {
    return $x;
  } else {
    return $y;
  }
}

# ------------------------------------------------------------------------

use constant {
	      JUNK_VALUE => 0,
	      COAL_VALUE => 1,
	      BRONZE_VALUE => 2,
	      SILVER_VALUE => 3,
	      GOLD_VALUE => 4
	     };

use constant {
	      COLOUR_JUNK => "255 255 255",   #FFFFFF
	      COLOUR_COAL => "0 0 0",	      #000000
	      COLOUR_BRONZE => "205 127 50",  #CD7F32
	      COLOUR_SILVER => "192 192 192", #C0C0C0
	      COLOUR_GOLD => "255 215 0"      #FFD700
	     };

sub parse_grade_name {
  my ($ss) = @_;
  my $s = lc($ss);

  return COAL_VALUE if ($s eq "black");
  return COAL_VALUE if ($s eq "coal");
  return BRONZE_VALUE if ($s eq "bronze");
  return SILVER_VALUE if ($s eq "silver");
  return GOLD_VALUE if ($s eq "gold");
  return GOLD_VALUE if ($s eq "golden");

  die "Unknown grade name: <<$ss>>,";
}

sub unparse_grade_name {
  my ($grade) = @_;

  return "coal" if ($grade == COAL_VALUE);
  return "bronze" if ($grade == BRONZE_VALUE);
  return "silver" if ($grade == SILVER_VALUE);
  return "gold" if ($grade == GOLD_VALUE);

  die "Unknown grade: <<$grade>>,";
}

sub colorize_grade {
  my ($grade) = @_;

  return COLOUR_JUNK if ($grade == JUNK_VALUE);
  return COLOUR_COAL if ($grade == COAL_VALUE);
  return COLOUR_BRONZE if ($grade == BRONZE_VALUE);
  return COLOUR_SILVER if ($grade == SILVER_VALUE);
  return COLOUR_GOLD if ($grade == GOLD_VALUE);

  die "Unknown grade: <<$grade>>,";
}

# ------------------------------------------------------------------------
# Process the command line
# ------------------------------------------------------------------------

use File::Basename;
use Getopt::Std;

our $opt_h;
our $opt_o;
our $opt_r;
our $opt_s = "silver";
our $opt_x = "silver";

my $progname = basename($0);

sub usage {

  print STDERR "Usage: $progname [options] simple-pg.gff tag1:annotation1.gff tag2:annotation2.gff ...\n";

  print STDERR "-h - print help\n";
  print STDERR "-o OUTPUT.gff - write results to OUTPUT.gff\n";
  print STDERR "-r OUTPUT.txt - write report to OUTPUT.txt\n";
  print STDERR "-s GRADE - Minimum ORF grade required for CDS support [$opt_s]\n";
  print STDERR "-x GRADE - Mininum peptide grade required for CDS extension [$opt_x]\n";

  exit(@_);
}

my $stat = getopts('ho:r:s:x:');
if (!$stat) {
  usage(1);
}
if ($opt_h) {
  usage();
}

if ( scalar(@ARGV) < 2 ) {
  usage(1);
}

$opt_s = parse_grade_name($opt_s);
$opt_x = parse_grade_name($opt_x);


# ------------------------------------------------------------------------

my $out_fh;
if ( !defined($opt_o) || $opt_o eq "-" ) {
  $out_fh = \*STDOUT;
} else {
  open($out_fh,">",$opt_o) || die "Cannot open for writing: <<$opt_o>>,";
}
my $report_fh;
if ($opt_r) {
  open($report_fh,">",$opt_r) || die "Cannot open for writing: <<$opt_r>>,\n";
}

sub report {
  print STDERR @_, "\n";
  if ($report_fh) {
    print $report_fh @_, "\n";
  }
}

# ------------------------------------------------------------------------

sub stringify_count {
  my ($n,$tag,$total) = @_;
  my $s = "";
  if ($tag ne "") {
    $s .= $tag.": ";
  }
  if (defined($total) && $total > 0) {
    $s .= sprintf("%d (%.2f%%)",
		  $n, ((100 * $n) / $total));
  } else {
    $s .= sprintf("%d", $n);
  }
  return $s;
}

sub stringify_counts {
  my ($grade_key,$a,$total) = @_;

  my $s1 = stringify_count(scalar(@$a),"",$total);

  if (!defined($total)) {
    $total = scalar(@$a);
  }
  my @counts = ( 0,		#junk
		 0,		#coal
		 0,		#bronze
		 0,		#silver
		 0);		#gold
  foreach my $entry (@$a) {
    my $grade = $entry->{attributes}->{$grade_key};
    $counts[$grade]++;
  }

  my @tags = ( "junk", "coal", "bronze", "silver", "gold" );
  my @l;
  for (my $i = scalar(@counts)-1; $i>=0; $i--) {
    if ($counts[$i] == 0) {
      next;
    }
    push @l, stringify_count($counts[$i], $tags[$i], $total);
  }
  my $s2 = "";
  if (scalar(@l) > 0) {
    $s2 = " (".join(", ", @l).")";
  }
  return $s1.$s2;
}


# ------------------------------------------------------------------------
# Read input gff files
# ------------------------------------------------------------------------

use constant { NO_VALUE => ";no-value;" };

sub parse_gff_attributes {
  my ($raw_attributes) = @_;
  my $attributes = {};
  foreach my $key_val (split(/; */,$raw_attributes)) {
    my ($key,$val);
    if ( $key_val =~ /^([^=]+)=(.*)/ ) {
      ($key,$val) = ($1,$2);
    } else {
      ($key,$val) = ($key_val, NO_VALUE);
    }
    $attributes->{$key} = $val;
  }
  return $attributes;
}

sub read_gff {
  my ($filename) = @_;

  report("## reading: $filename");

  my @entries;
  open(my $fh, "<", $filename) || die "Cannot open <<$filename>>,";
  while (<$fh>) {
    if ( /^#/ ) {
      next;
    }
    chomp;
    my ($accession,$source,$feature,$start,$end,
	$score,$strand,$frame,$raw_attributes) = split(/\t/,$_);
    my $attributes = parse_gff_attributes($raw_attributes);

    $accession =~ s/\.[0-9]+$//;

    (defined($source)) || die;
    my $entry = {
		 accession => $accession,
		 source => $source,
		 feature => $feature,
		 start => $start,
		 end => $end,
		 score => $score,
		 strand => $strand,
		 frame => $frame,
		 attributes => $attributes,
		 accession => $accession,
		};

    push @entries, $entry;
  }
  close $fh;

  report(sprintf("## ... %d features",scalar(@entries)));

  return @entries;
}

# ------------------------------------------------------------------------
# Filtering proteomics features
# ------------------------------------------------------------------------
report("#\n# Filtering proteomics features");

my $proteomics_gff = shift @ARGV;
my @proteomics_entries = read_gff($proteomics_gff);


## Find the supported ORF's
my %supported_orfs;
foreach my $entry ( @proteomics_entries ) {
  # Only ORFS
  if ( $entry->{feature} ne "ORF" ) {
    next;
  }
  my $grade = parse_grade_name($entry->{attributes}->{grade});
  if ( $opt_s > $grade ) {
    next;
  }
  my $name = $entry->{attributes}->{name};
  (defined($name)) || die;
  (!defined($supported_orfs{$name})) || die;
  my $o = { %$entry } ;
  $o->{attributes} = {grade => $grade};
  $o->{__peptides__} = [];
  $supported_orfs{$name} = $o;
}

my @supported_orfs = values(%supported_orfs);

report("## Supported ORF's: ",
       stringify_counts("grade",[@supported_orfs],scalar(@supported_orfs)));

## Find the supporting peptides;
foreach my $entry ( @proteomics_entries ) {
  # looking for plain_peptide's
  if ( $entry->{feature} ne "plain_peptide" ) {
    next;
  }
  # Only consider peptides that are good enough to support extension
  my $grade = parse_grade_name($entry->{attributes}->{grade});
  if ( $opt_x > $grade ) {
    next;
  }
  # Get the matching ORF
  my $name = $entry->{attributes}->{orf_name};
  my $orf = $supported_orfs{$name};
  if (!$orf) {
    # The matching ORF didn't make the grade...
    next;
  }

  my $o = { %$entry } ;
  $o->{attributes} = {grade => $grade};
  push @{$orf->{__peptides__}}, $o;
}

# ------------------------------------------------------------------------
# Routines to reconcile support ORF's with annotations
# ------------------------------------------------------------------------

sub accession_index {
  my ($spatial_indices,$accession) = @_;
  my $spatial_index = $spatial_indices->{$accession};
  if (!$spatial_index) {
    $spatial_index = $spatial_indices->{$accession} =
      Set::IntervalTree->new;
  }
  return $spatial_index;
}

sub frame_index {
  my ($e) = @_;
  my $strand = $e->{strand};
  my $start = $e->{start};
  my $end = $e->{end};
  if ($strand eq "+") {
    return  $start%3+1;
  } else {
    return -($end%3+1);
  }
}

sub overlapping_entries {
  my ($target,@queries) = @_;
  my $t_start = $target->{start};
  my $t_end = $target->{end};
  # my $t_frame_index = frame_index($target);
  my @results;
  foreach my $q_entry ( @queries ) {
    # ($t_frame_index == frame_index($q_entry)) || die;
    my $q_start = $q_entry->{start};
    my $q_end = $q_entry->{end};
    if ( ( $t_start <= $q_start && $q_start <= $t_end )
	 || ( $t_start <= $q_end && $q_end <= $t_end ) ) {
      push @results, $q_entry;
    }
  }
  return @results;
}

sub reconcile {
  my ($db_annotation_gff) = @_;
  report("#\n# Reconciling $db_annotation_gff");

  # ------------------------------------------------------------------------
  # Index annotation entries
  # ------------------------------------------------------------------------

  report("## Indexing annotation entries");

  my @l = split(/:/,$db_annotation_gff);
  my ($db_prefix,$annotation_gff);
  if (scalar(@l) == 1 ) {
    ($db_prefix,$annotation_gff) = ("",$l[0]);
  } elsif (scalar(@l) == 2) {
    ($db_prefix,$annotation_gff) = ($l[0]."_",$l[1]);
  } else {
    die "Invalid argument <<$db_annotation_gff>>,";
  }

  my $_annotation_edit = $db_prefix."annotation_edit";
  my $_supported_cds = $db_prefix."supported_CDS";
  my $_locus_tag = $db_prefix."locus_tag";
  my $_check_manually = $db_prefix."check_manually";
  my $_coding_pseudogene = $db_prefix."coding_pseudogene";
  my $_multi_cds = $db_prefix."multi_CDS";
  my $_novel_cds = $db_prefix."novel_CDS";
  my $_short_cterm = $db_prefix."short_Cterm";
  my $_short_nterm = $db_prefix."short_Nterm";


  my @annotation_entries = read_gff($annotation_gff);

  # store and retrieve half closed intervals [lb,ub)
  my $spatial_indices = {};

  my $num_CDSs = 0;
  foreach my $entry (@annotation_entries) {
    my $accession = $entry->{accession};
    my $start = $entry->{start};
    my $end = $entry->{end};
    ($start < $end) || die "<<$_>>,";
    accession_index($spatial_indices,$accession)->insert($entry,$start,$end+1);
    if ($entry->{feature} eq "CDS") {
      $num_CDSs++;
    }
  }

  # ------------------------------------------------------------------------
  # Reconciling supported orfs with annotations
  # ------------------------------------------------------------------------
  report("## Reconciling supported orfs with annotations");

  my @supported_CDSs = ();
  my @novel_CDSs = ();
  my @coding_pseudogenes = ();
  my @multi_CDS_ORFs = ();
  my @short_nterminal_CDSs = ();
  my @short_cterminal_CDSs = ();

  foreach my $orf_entry ( @supported_orfs ) {

    my $orf_accession = $orf_entry->{accession};
    my $orf_start = $orf_entry->{start};
    my $orf_end = $orf_entry->{end};
    my $orf_strand = $orf_entry->{strand};
    my $orf_frame_index = frame_index($orf_entry);
    my $orf_evidence = $orf_entry->{attributes}->{grade};
    my @orf_peptides = @{$orf_entry->{__peptides__}};

    # --------------------------------------------------
    # Find and note the n-terminal position of the
    # n-terminal-most peptide
    # --------------------------------------------------

    my ($min_start,$max_end);
    foreach my $peptide_entry ( @orf_peptides ) {
      my $peptide_start = $peptide_entry->{start};
      my $peptide_end = $peptide_entry->{end};
      if (!defined($min_start) || $peptide_start < $min_start) {
	$min_start = $peptide_start;
      }
      if (!defined($max_end) || $peptide_end > $max_end) {
	$max_end = $peptide_end;
      }
    }
    defined($min_start) || die;
    defined($max_end) || die;
    if ($orf_strand eq "+") {
      $orf_entry->{attributes}->{observed_Nterm} = $min_start;
    } else {
      $orf_entry->{attributes}->{observed_Nterm} = $max_end;
    }

    # --------------------------------------------------
    # find entries in frame
    # --------------------------------------------------

    my $overlapping_annotations =
      accession_index($spatial_indices,$orf_accession)->fetch($orf_start,$orf_end+1);
    (defined($overlapping_annotations)) || die;
    my @in_frame_annotations;
    foreach my $annot_entry (@{$overlapping_annotations}) {
      if ( $orf_frame_index != frame_index($annot_entry) ) {
	next;
      }
      push @in_frame_annotations, $annot_entry;
    }

    # --------------------------------------------------
    # coding_pseuogene == gene or CDS with /pseudo
    # --------------------------------------------------

    my $found_pseudo = FALSE;
    foreach my $annot_entry (@in_frame_annotations) {
      if ( $annot_entry->{feature} ne "gene" && $annot_entry->{feature} ne "CDS" ) {
	next;
      }
      if (!defined($annot_entry->{attributes}->{pseudo})) {
	next;
      }
      $found_pseudo = TRUE;
    }
    if ($found_pseudo) {
      $orf_entry->{attributes}->{$_coding_pseudogene} = $orf_evidence;
      push @coding_pseudogenes, $orf_entry;
      $orf_entry->{attributes}->{$_annotation_edit} = NO_VALUE;
    }


    # --------------------------------------------------
    # find the CDS annotation(s)
    # --------------------------------------------------

    my @CDS_annotations;
    foreach my $annot_entry (@in_frame_annotations) {
      if ( $annot_entry->{feature} ne "CDS" ) {
	next;
      }
      push @CDS_annotations, $annot_entry;
    }

    # --------------------------------------------------
    # if no CDS's, then mark ORF as novel
    # --------------------------------------------------

    if ( scalar(@CDS_annotations) == 0 ) {
      $orf_entry->{attributes}->{$_novel_cds} = $orf_evidence;
      push @novel_CDSs, $orf_entry;
      $orf_entry->{attributes}->{$_annotation_edit} = NO_VALUE;
      next;
    }

    # --------------------------------------------------
    # Multi-CDS ORFs == contains >1 supported CDS (e.g. P.ananatis LMG20103)
    # --------------------------------------------------

    if ( scalar(@CDS_annotations) > 1 ) {
      $orf_entry->{attributes}->{$_multi_cds} = $orf_evidence;
      push @multi_CDS_ORFs, $orf_entry;
      $orf_entry->{attributes}->{$_check_manually} = NO_VALUE;
      next;
    }


    # --------------------------------------------------
    # single CDS ORF's
    # --------------------------------------------------

    (scalar(@CDS_annotations) == 1) || die;
    $orf_entry->{attributes}->{$_supported_cds} = $orf_evidence;
    push @supported_CDSs, $orf_entry;


    my ($annot_entry) = @CDS_annotations;

    my $locus_tag = $annot_entry->{attributes}->{locus_tag};
    (defined($locus_tag)) || die;
    $orf_entry->{attributes}->{$_locus_tag} = $locus_tag;

    # --------------------------------------------------
    # CDS's that need lengthening
    # --------------------------------------------------

    # make skeletal entries for the regions within the ORF "left" and
    # "right" of the existing CDS annotation.
    my $left_entry = {
		      start => $orf_entry->{start},
		      end => $annot_entry->{start}-1,
		      strand => $orf_entry->{strand}
		     };
    my $right_entry = {
		       start => $annot_entry->{end}+1,
		       end => $orf_entry->{end},
		       strand => $orf_entry->{strand}
		      };

    # map "left" and "right" to "N-terminal" and "C-terminal" based on the strand
    my ($nterm_entry,$cterm_entry);
    if ($orf_entry->{strand} eq "+") {
      ($nterm_entry,$cterm_entry) = ($left_entry,$right_entry);
    } else {
      ($orf_entry->{strand} eq "-") || die;
      ($nterm_entry,$cterm_entry) = ($right_entry,$left_entry);
    }

    # Find "supported" peptides that overlap the non-empty n-terminal region
    if ($nterm_entry->{start} <= $nterm_entry->{end}) {
      my @l = overlapping_entries($nterm_entry,@orf_peptides);
      if (scalar(@l) > 0) {
	my $peptide_evidence = JUNK_VALUE;
	foreach my $peptide_entry ( @l ) {
	  if ( $peptide_evidence < $peptide_entry->{attributes}->{grade} ) {
	    $peptide_evidence = $peptide_entry->{attributes}->{grade};
	  }
	}
	$orf_entry->{attributes}->{$_short_nterm} = $peptide_evidence;
	push @short_nterminal_CDSs, $orf_entry;
	$orf_entry->{attributes}->{$_annotation_edit} = NO_VALUE;
      }
    }

    # Find "supported" peptides that overlap the non-empty c-terminal region
    if ($cterm_entry->{start} <= $cterm_entry->{end}) {
      my @l = overlapping_entries($cterm_entry,@orf_peptides);
      if (scalar(@l) > 0) {
	my $peptide_evidence = JUNK_VALUE;
	foreach my $peptide_entry ( @l ) {
	  if ( $peptide_evidence < $peptide_entry->{attributes}->{grade} ) {
	    $peptide_evidence = $peptide_entry->{attributes}->{grade};
	  }
	}
	$orf_entry->{attributes}->{$_short_cterm} = $peptide_evidence;
	push @short_cterminal_CDSs, $orf_entry;
	$orf_entry->{attributes}->{$_annotation_edit} = NO_VALUE;
      }
    }

  }

  # ------------------------------------------------------------------------
  # Output summary.
  # ------------------------------------------------------------------------

  report("## $_supported_cds\'s (of $num_CDSs): ",
	 stringify_counts($_supported_cds,[@supported_CDSs], $num_CDSs));
  report("## $_short_nterm\'s: ",
	 stringify_counts($_short_nterm,[@short_nterminal_CDSs]));
  report("## $_novel_cds\'s: ",
	 stringify_counts($_novel_cds,[@novel_CDSs]));
  report("## $_coding_pseudogene\'s: ",
	 stringify_counts($_coding_pseudogene,[@coding_pseudogenes]));
  if (scalar(@multi_CDS_ORFs) > 0) {
    report("## >>> $_multi_cds\'s <<<: ",
	   stringify_counts($_multi_cds,[@multi_CDS_ORFs]));
  }
  if (scalar(@short_cterminal_CDSs) > 0) {
    report("## >>> $_short_cterm\'s <<<: ",
	   stringify_counts($_short_cterm,[@short_cterminal_CDSs]));
  }

  foreach my $e ( @supported_orfs ) {
    my $attrs = $e->{attributes};
    foreach my $k (

		   $_annotation_edit,
		   $_check_manually,
		   $_coding_pseudogene,
		   $_multi_cds,
		   $_novel_cds,
		   $_supported_cds,

		  ) {
      if (defined($attrs->{$k})) {
	$attrs->{$k} = NO_VALUE;
      }
    }

    foreach my $k (
		   $_short_cterm,
		   $_short_nterm,
		  ) {
      if (defined($attrs->{$k})) {
	$attrs->{$k} = unparse_grade_name($attrs->{$k});
      }
    }
  }

}

# ------------------------------------------------------------------------
# Perform reconciliations
# ------------------------------------------------------------------------

foreach my $db_annotation_gff ( @ARGV ) {
  reconcile($db_annotation_gff);
}

# ------------------------------------------------------------------------
# Output results.
# ------------------------------------------------------------------------

report("#\n# Output results");

# ------------------------------------------------------------------------

sub unparse_gff_attributes {
  my ($attributes) = @_;
  my @l;
  foreach my $key (sort (keys %$attributes) ) {
    my $val = $attributes->{$key};
    if (!defined($val)) {
      next;
    } elsif ( $val eq NO_VALUE ) {
      push @l, $key;
    } else {
      push @l, $key."=".$val;
    }
  }
  return join(";",@l);
}

print $out_fh "##gff-version 3\n";

my @output_entries = @supported_orfs;
# determinist output rocks.
@output_entries = sort { $a->{accession} cmp $b->{accession}
			   || $a->{start} <=> $b->{start}
			   || $a->{end} <=> $b->{end}
			   || $a->{strand} cmp $b->{strand}
			 } @output_entries;

foreach my $e ( @output_entries ) {
  my $attrs = $e->{attributes};
  my $grade = $attrs->{grade};
  $attrs->{grade} = unparse_grade_name($grade);
  $attrs->{colour} = colorize_grade($grade);

  print $out_fh
    join("\t",
	 $e->{accession},
	 $e->{source},
	 $e->{feature},
	 $e->{start},
	 $e->{end},
	 $e->{score},
	 $e->{strand},
	 $e->{frame},
	 unparse_gff_attributes($e->{attributes})
	),"\n";
}

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------
report("#\n# Done");

if ( defined($opt_o) && $opt_o ne "-" ) {
  close $out_fh;
}
if ( $opt_r ) {
  close($report_fh);
}
