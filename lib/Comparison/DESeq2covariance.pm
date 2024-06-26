#!/usr/bin/perl
package Comparison::DESeq2covariance;

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

  copy(dirname(__FILE__) . "/../CQS/countTableVisFunctions.R",  "$result_dir/countTableVisFunctions.R");

  #print("writeParameterSampleFile in DESeq2\n");
  writeParameterSampleFile( $config, $section, $result_dir, 1 );

  my $comparisons = get_raw_files( $config, $section );
  my @comparison_names = keys %{$comparisons};
  
  my $groups = get_raw_files( $config, $section, "groups" );

  my $covariate_file = get_option_file($config, $section, "covariance_file", 1);
  my ($cov_table, $cov_names) = read_table($covariate_file, get_option($config, $section, "covariance_name_index", 0));

  #print(Dumper($cov_table));

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
  if ( is_hash($countfiles) ) {
    if ( scalar( keys %$countfiles ) > 1 ) {
      $single_count_file = 0;
    }
  }

  my $designfilename = "${task_name}.define";
  my $designfile     = "$result_dir/$designfilename";
  open( my $df, ">$designfile" ) or die "Cannot create $designfile";
  print $df "ComparisonName\tCountFile\tConditionFile\tReferenceGroupName\tSampleGroupName\tComparisonTitle\tdesignFormula\tcontrast\tcollapse_by\tpairOnlyCovariant\n";

  for my $comparisonIndex ( 0 .. $#comparison_names ) {
    my $comparison_name = $comparison_names[$comparisonIndex];
    my $comparisonTitle = $comparisonTitles->[$comparisonIndex];
    $first++;

    my $comp_def = $comparisons->{$comparison_name};
    print "$comparison_name\n";

    my $group_names;
    my $covariate_names;
    my $contrast="";
    my $designFormula="";
    my $collapse_by="";
    my $pairOnlyCovariant="";
    if (is_array($comp_def)){
      $group_names = $comp_def;
      $covariate_names = [];
    }elsif (is_hash($comp_def)){
      $group_names = $comp_def->{"groups"};
      $covariate_names = $comp_def->{covariances};
      if (defined $comp_def->{"designFormula"}) {
        my $formula = $comp_def->{"designFormula"};
        if(is_array($formula)){
          $designFormula = $formula->[0];
        }else{
          $designFormula = $formula;
        }
      }
      if (defined $comp_def->{"contrast"}) {
        my $cons = $comp_def->{"contrast"};
        if(is_array($cons)){
          $contrast = $cons->[0];
        }else{
          $contrast = $cons;
        }
      }
      if (defined $comp_def->{"collapse_by"}) {
        $collapse_by = $comp_def->{collapse_by};
      }
      if (defined $comp_def->{"pairOnlyCovariant"}) {
        $pairOnlyCovariant = $comp_def->{pairOnlyCovariant};
      }
    }else{
      die("Comparison $comparison_name is not defined correctly, check your configuration");
    }

    if(!is_array($group_names)){
      die("groups of comparison $comparison_name is not array, check your configuration");
    }

    if ( scalar(@$group_names) != 2 ) {
      die "Comparison of $comparison_name should contains and only contains two groups!";
    }

    for my $covariance (@$covariate_names){
      if (not defined $cov_names->{$covariance}){
        die "Cannot find covariate $covariance of comparison $comparison_name in covarite file $covariate_file";
      }
    }
    my @covariate_keys = @$covariate_names;

    #print( Dumper(@group_names) );

    my $g1 = $group_names->[0];
    my $g2 = $group_names->[1];
    die "cannot find group $g1 " if !defined( $groups->{$g1} );
    die "cannot find group $g2 " if !defined( $groups->{$g2} );
    my @s1 = sort @{ $groups->{$g1} };
    my @s2 = sort @{ $groups->{$g2} };
    
    for my $s11 (@s1){
      if (not defined $cov_table->{$s11}){
        die "Cannot find sample $s11 of group $g1 of comparison $comparison_name in covariance file $covariate_file";
      }
    }
    
    for my $s22 (@s2){
      if (not defined $cov_table->{$s22}){
        die "Cannot find sample $s22 of group $g2 of comparison $comparison_name in covariance file $covariate_file";
      }
    }

    my $total_sample_count = scalar(@s1) + scalar(@s2);

    my $filename = "${comparison_name}.design";

    my $cdfile = $result_dir . "/$filename";
    open( my $cd, ">$cdfile" ) or die "Cannot create $cdfile";
    if ( scalar(@covariate_keys) > 0 ) {
      print $cd "Sample\tCondition\t", join( "\t", @covariate_keys ), "\n";
    }
    else {
      print $cd "Sample\tCondition\n";
    }
    for my $i ( 0 .. $#s1 ) {
      my $sname = $s1[$i];
      print $cd "${sname}\t${g1}";
      if ( scalar(@covariate_keys) > 0 ) {
        for my $key (@covariate_keys) {
          print $cd "\t" . $cov_table->{$sname}{$key};
        }
      }
      print $cd "\n";
    }
    for my $i ( 0 .. $#s2 ) {
      my $sname = $s2[$i];
      print $cd "${sname}\t${g2}";
      if ( scalar(@covariate_keys) > 0 ) {
        for my $key (@covariate_keys) {
          print $cd "\t" . $cov_table->{$sname}{$key};
        }
      }
      print $cd "\n";
    }
    close $cd;

    my $curcountfile = $single_count_file ? $countfile : $countfiles->{$comparison_name};
    if ( !defined $curcountfile ) {
      print Dumper($countfiles);
      die "Count file not found for comparison of $comparison_name!";
    }
    if ( ref $curcountfile eq ref [] ) {
      $curcountfile = $curcountfile->[0];
    }
    print $df "$comparison_name\t$curcountfile\t$cdfile\t$g1\t$g2\t$comparisonTitle\t$designFormula\t$contrast\t$collapse_by\t$pairOnlyCovariant\n";
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
