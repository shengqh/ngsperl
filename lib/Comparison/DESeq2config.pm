#!/usr/bin/perl
package Comparison::DESeq2config;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use File::Copy;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use Comparison::DESeq2;

our @ISA = qw(Comparison::DESeq2);

my $directory;

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_de2c";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $comparisons = get_raw_files( $config, $section );
  my @comparison_names = keys %{$comparisons};

  my $comparisonAttributes = get_raw_files_attributes( $config, $section );
  my $comparisonTitles = \@comparison_names;
  if ( defined $comparisonAttributes->{".order"} && defined $comparisonAttributes->{".col"} ) {
    $comparisonTitles = $comparisonAttributes->{".col"};
  }

  my $totalPair = scalar(@comparison_names);
  if ( 0 == $totalPair ) {
    die "No pair defined!";
  }

  my $countfile = parse_param_file( $config, $section, "countfile", 1 );
  my $rtemplate = dirname(__FILE__) . "/DESeq2.r";
  if ( !-e $rtemplate ) {
    die "File not found : " . $rtemplate;
  }

  my $showLabelInPCA    = get_option( $config, $section, "show_label_PCA",       1 );
  my $showDEGeneCluster = get_option( $config, $section, "show_DE_gene_cluster", 0 );
  my $pvalue            = get_option( $config, $section, "pvalue",               0.05 );
  my $foldChange        = get_option( $config, $section, "fold_change",          2.0 );
  my $minMedianInGroup  = get_option( $config, $section, "min_median_read",      0 );
  my $addCountOne       = get_option( $config, $section, "add_count_one",        0 );
  my $usePearsonInHCA   = get_option( $config, $section, "use_pearson_in_hca",   0 );

  my $top25only                 = get_option( $config, $section, "top25only",                    0 );
  my $detectedInBothGroup       = get_option( $config, $section, "detected_in_both_group",       0 );
  my $performWilcox             = get_option( $config, $section, "perform_wilcox",               0 );
  my $useRawPvalue              = get_option( $config, $section, "use_raw_p_value",              0 );
  my $textSize                  = get_option( $config, $section, "text_size",                    11 );
  my $transformTable            = get_option( $config, $section, "transform_table",              0 );
  my $exportSignificantGeneName = get_option( $config, $section, "export_significant_gene_name", 0 );
  my $cooksCutoff               = get_option( $config, $section, "cooksCutoff",                  'DEFAULT' );
  my $independentFiltering      = get_option( $config, $section, "independentFiltering", 'TRUE' );
  
  my $rCode = get_option( $config, $section, "rCode", "" );

  my $libraryFile = parse_param_file( $config, $section, "library_file", 0 );
  my $libraryKey;
  if ( defined $libraryFile ) {
    $libraryKey = get_option( $config, $section, "library_key", "TotalReads" );
  }

  my $suffix = $self->getSuffix( $top25only, $detectedInBothGroup, $minMedianInGroup, $useRawPvalue, $pvalue );

  my $first = 0;

  my $countfiles = defined $config->{$section}{"countfile"} ? {} : get_raw_files( $config, $section, "countfile" );
  my $single_count_file = 1;
  if ( ref $countfiles eq ref {} ) {
    if ( scalar( keys %$countfiles ) > 1 ) {
      $single_count_file = 0;
    }
  }

  my $designfilename = "${task_name}.define";
  my $designfile     = "$result_dir/$designfilename";
  open( my $df, ">$designfile" ) or die "Cannot create $designfile";
  print $df "ComparisonName\tCountFile\tConditionFile\tReferenceGroupName\tSampleGroupName\tComparisonTitle\n";

  for my $comparisonIndex ( 0 .. $#comparison_names ) {
    my $comparison_name = $comparison_names[$comparisonIndex];
    my $comparisonTitle = $comparisonTitles->[$comparisonIndex];
    $first++;

    my $pairOnlyCovariant="";

    my $design_file = $comparisons->{$comparison_name};

    my $filename = "${comparison_name}.design";
    my $cdfile = $result_dir . "/$filename";

    copy($design_file,$cdfile) or die "Copy failed: $!";

    my $curcountfile = $single_count_file ? $countfile : $countfiles->{$comparison_name};
    if ( !defined $curcountfile ) {
      print Dumper($countfiles);
      die "Count file not found for comparison of $comparison_name!";
    }
    if ( ref $curcountfile eq ref [] ) {
      $curcountfile = $curcountfile->[0];
    }
    print $df "$comparison_name\t$curcountfile\t$cdfile\tcontrol\tcase\t$comparisonTitle\n";
  }
  close($df);

  my $rfile = $result_dir . "/${task_name}.r";
  open( my $rf, ">$rfile" )     or die "Cannot create $rfile";
  open( my $rt, "<$rtemplate" ) or die $!;
  print $rf "
rootdir<-\"$result_dir\"
inputfile<-\"$designfilename\" 

pvalue<-$pvalue
useRawPvalue<-$useRawPvalue
foldChange<-$foldChange
minMedianInGroup<-$minMedianInGroup
  
detectedInBothGroup<-$detectedInBothGroup
showLabelInPCA<-$showLabelInPCA
showDEGeneCluster<-$showDEGeneCluster
addCountOne<-$addCountOne
usePearsonInHCA<-$usePearsonInHCA
top25only<-$top25only
performWilcox<-$performWilcox
textSize<-$textSize
transformTable<-$transformTable
exportSignificantGeneName<-$exportSignificantGeneName
thread<-$thread

independentFiltering<-$independentFiltering

$rCode
";

  if ( $cooksCutoff ne "DEFAULT" ) {
    print $rf "cooksCutoff<-$cooksCutoff
";
  }

  if ( defined $libraryFile ) {
    print $rf "
libraryFile<-\"$libraryFile\"
libraryKey<-\"$libraryKey\"
";
  }
  
  print $rf "#predefined_condition_end\n";

  while (<$rt>) {
    if ( $_ !~ 'predefined_condition_end' ) {
      next;
    }
    last;
  }
  while (<$rt>) {
    s/\r|\n//g;
    print $rf $_, "\n";
  }
  close $rt;
  close $rf;

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $lastComparisonName = $comparison_names[-1];
  my $final_file         = $lastComparisonName . $suffix . "_DESeq2_volcanoPlot.png";
  my $pbs                = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  print $pbs "R --vanilla -f $rfile \n";

  $self->close_pbs( $pbs, $pbs_file );
}

1;
