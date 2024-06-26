#!/usr/bin/perl
package Comparison::DESeq2contrast;

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
use CQS::UniqueTask;
use CQS::StringUtils;

our @ISA = qw(CQS::UniqueTask);

my $directory;

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_de2";
  bless $self, $class;
  return $self;
}

sub getSuffix {
  my ( $self, $top25only, $detectedInBothGroup, $minMedianInGroup, $useRawPvalue, $pvalue ) = @_;
  my $suffix = "";
  if ($top25only) {
    $suffix = $suffix . "_top25";
  }
  if ($detectedInBothGroup) {
    $suffix = $suffix . "_detectedInBothGroup";
  }
  if ( $minMedianInGroup > 0 ) {
    $suffix = $suffix . "_min${minMedianInGroup}";
  }
  if ($useRawPvalue) {
    $suffix = $suffix . "_pvalue" . $pvalue;
  }
  else {
    $suffix = $suffix . "_fdr" . $pvalue;
  }
  return $suffix;
}

sub output_report {
  my ( $task_name, $pbs_dir, $result_dir ) = @_;

  copy(dirname(__FILE__) . "/DESeq2.rmd",  "$result_dir/DESeq2.rmd");
  copy(dirname(__FILE__) . "/../CQS/reportFunctions.R",  "$result_dir/reportFunctions.R");

  my $rmd_command = "Rscript --vanilla -e \"library('rmarkdown');rmarkdown::render('DESeq2.rmd',output_file='${task_name}.deseq2.html')\"";
  my $report_sh = "$pbs_dir/report.sh";
  open( my $rs, ">$report_sh" ) or die $!;
  print $rs "
cd $result_dir

$rmd_command
";
  close($rs);

  return($rmd_command);
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  copy(dirname(__FILE__) . "/../CQS/countTableVisFunctions.R",  "$result_dir/countTableVisFunctions.R");

  #print("writeParameterSampleFile in DESeq2\n");
  writeParameterSampleFile( $config, $section, $result_dir, 1 );

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

  my $groups = get_raw_files( $config, $section, "groups" );

  my $countfile = parse_param_file( $config, $section, "countfile", 1 );
  my $rtemplate = dirname(__FILE__) . "/DESeq2contrast.r";
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

  writeParameterSampleFile( $config, $section, $result_dir, 1 );
  
  my $rCode = get_option( $config, $section, "rCode", "" );

  my $libraryFile = parse_param_file( $config, $section, "library_file", 0 );
  my $libraryKey;
  if ( defined $libraryFile ) {
    $libraryKey = get_option( $config, $section, "library_key", "TotalReads" );
  }

  my $suffix = $self->getSuffix( $top25only, $detectedInBothGroup, $minMedianInGroup, $useRawPvalue, $pvalue );

  my %tpgroups = ();
  for my $group_name ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$group_name} };
    $tpgroups{$group_name} = "\"" . join( "\",\"", @samples ) . "\"";
  }

  my $first = 0;

  my $countfiles = defined $config->{$section}{"countfile"} ? {} : get_raw_files( $config, $section, "countfile" );
  my $single_count_file = 1;
  if ( ref $countfiles eq ref {} ) {
    if ( scalar( keys %$countfiles ) > 1 ) {
      $single_count_file = 0;
    }
  }

  my $cov_table = undef;
  my $cov_names = undef;
  my $covariate_file = $config->{$section}{"covariance_file"};
  if(defined $covariate_file){
    ($cov_table, $cov_names) = read_table($covariate_file, get_option($config, $section, "covariance_name_index", 0));
  }

  my $designfilename = "${task_name}.define";
  my $designfile     = "$result_dir/$designfilename";
  open( my $df, ">$designfile" ) or die "Cannot create $designfile";
  print $df "ComparisonName\tCountFile\tConditionFile\tReferenceGroupName\tSampleGroupName\tComparisonTitle\tpairOnlyCovariant\tdesignFormula\tcontrast\n";

  for my $comparisonIndex ( 0 .. $#comparison_names ) {
    my $comparison_name = $comparison_names[$comparisonIndex];
    my $comparisonTitle = $comparisonTitles->[$comparisonIndex];
    $first++;

    my $covariances = {};
    my $pairOnlyCovariant="";

    my $gNames = $comparisons->{$comparison_name};
    if ( !( ref $gNames eq ref {} ) ) {
      die "Definition of " . $comparison_name . " should be hash table to include groups and contrast!";
    }

    my @group_names = @{ $gNames->{groups} };
    my $contrast = $gNames->{contrast};
    if ( !defined $contrast ) {
      die "Definition of " . $comparison_name . " should have contrast defined!";
    }
    if(is_array($contrast)) {
      die "contrast should be string, current is array for comparison " . $comparison_name;
    }

    my $designFormula = "";
    if(defined $gNames->{"designFormula"}){
      $designFormula = $gNames->{"designFormula"};
      if(is_array($designFormula)){
        die "designFormula should be string, current is array for comparison " . $comparison_name;
      }
    }

    for my $key ( sort keys %$gNames ) {
      next if ( $key eq "groups" );
      next if ( $key eq "contrast" );
      next if ( $key eq "pairOnlyCovariant" );
      next if ( $key eq "designFormula" );
      next if ( $key eq "covariances" );
      
      $covariances->{$key} = $gNames->{$key};
    }

    my @all_samples = ();
    for my $gname (@group_names){
      if(!defined($groups->{$gname})){
        die "Cannot find group $gname in groups definition!";
      }
      my @samples = @{ $groups->{$gname} };
      push(@all_samples, @samples);
    }

    if(defined $gNames->{covariances}){
      die "no covariance_file defined, but having covariances in comparison definition." if !defined $covariate_file;
      my $covariance_names = $gNames->{covariances};
      for my $cov_name (@$covariance_names){
        if(!defined $cov_table->{$all_samples[0]}{$cov_name}){
          die "Cannot find covariance column $cov_name in $covariate_file!";
        }
        my $cov_values = [];
        for my $sample (@all_samples){
          push(@$cov_values, $cov_table->{$sample}{$cov_name});
        }

        $covariances->{$cov_name} = $cov_values;
      }
    }

    if (defined $gNames->{"pairOnlyCovariant"}) {
      $pairOnlyCovariant = $gNames->{pairOnlyCovariant};
    }

    my @covariances_keys = sort keys %$covariances;

    #print( Dumper(@group_names) );

    my $total_sample_count = 0;
    for my $gname (@group_names){
      if(!defined($groups->{$gname})){
        die "Cannot find group $gname in groups definition!";
      }

      my @samples = sort @{ $groups->{$gname} };
      $total_sample_count = $total_sample_count + scalar(@samples);
    }

    for my $key ( keys %$covariances ) {
      my $values = $covariances->{$key};

      if ( !( ref $values eq ref [] ) ) {
        die "Covariances of " . $key . " should be array reference!";
      }
      elsif ( scalar(@$values) != $total_sample_count ) {
        die "Number of covariance value of " . $key . " shoud be $total_sample_count !";
      }
    }    

    my $filename = "${comparison_name}.design";

    my $cdfile = $result_dir . "/$filename";
    open( my $cd, ">$cdfile" ) or die "Cannot create $cdfile";
    if ( scalar(@covariances_keys) > 0 ) {
      print $cd "Sample\tCondition\t", join( "\t", @covariances_keys ), "\n";
    }
    else {
      print $cd "Sample\tCondition\n";
    }

    my $cov_index = 0;
    for my $gname (@group_names){
      if(!defined($groups->{$gname})){
        die "Cannot find group $gname in groups definition!";
      }
      my @samples = @{ $groups->{$gname} };

      for my $sname ( @samples ) {
        print $cd "${sname}\t${gname}";
        if ( scalar(@covariances_keys) > 0 ) {
          for my $key (@covariances_keys) {
            print $cd "\t" . $covariances->{$key}[$cov_index];
          }
        }
        print $cd "\n";
        $cov_index = $cov_index + 1;
      }
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
    my $gcontrol = $group_names[ 0 ];
    #print($gcontrol, "\n");
    my @gsamples = @group_names[ 1 .. $#group_names ];
    #print(@gsamples, "\n");
    my $gsamples_str = join(",", @gsamples);

    print $df "$comparison_name\t$curcountfile\t$cdfile\t$gcontrol\t$gsamples_str\t$comparisonTitle\t$pairOnlyCovariant\t$designFormula\t$contrast\n";
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

  my $rmd_cmd = output_report($task_name, $pbs_dir, $result_dir);
  print $pbs "\n$rmd_cmd \n";

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern, $removeEmpty ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $comparisons = get_raw_files( $config, $section );
  my $minMedianInGroup          = get_option( $config, $section, "min_median_read",              0 );
  my $top25only                 = get_option( $config, $section, "top25only",                    0 );
  my $detectedInBothGroup       = get_option( $config, $section, "detected_in_both_group",       0 );
  my $performWilcox             = get_option( $config, $section, "perform_wilcox",               0 );
  my $exportSignificantGeneName = get_option( $config, $section, "export_significant_gene_name", 0 );
  my $useRawPvalue              = get_option( $config, $section, "use_raw_p_value",              0 );
  my $pvalue                    = get_option( $config, $section, "pvalue",                       0.05 );
  my $suffix = $self->getSuffix( $top25only, $detectedInBothGroup, $minMedianInGroup, $useRawPvalue, $pvalue );
  my $result = {};

  my $tasknameFiles = [$result_dir . "/${task_name}.define.DESeq2.version"];
  if ( scalar( keys %$comparisons ) > 1 ) {
    push @$tasknameFiles, $result_dir . "/${task_name}.define${suffix}_DESeq2_volcanoPlot.png";
  }
  my $filtered = filter_array( $tasknameFiles, $pattern, 1 );
  if ( scalar(@$filtered) > 0 || !$removeEmpty ) {
    $result->{$task_name} = $filtered;
  }

  for my $comparison_name ( sort keys %{$comparisons} ) {
    my @result_files = ();
    my $prefix       = $comparison_name . $suffix;

    my $final_file = $prefix . "_DESeq2_sig.csv";

    push( @result_files, $result_dir . "/${prefix}.csv" );
    push( @result_files, $result_dir . "/${prefix}_DESeq2.csv" );
    push( @result_files, $result_dir . "/${prefix}_DESeq2-vsd.csv" );
    push( @result_files, $result_dir . "/${prefix}_DESeq2_GSEA.rnk" );
    push( @result_files, $result_dir . "/${prefix}_DESeq2_sig.csv" );
    push( @result_files, $result_dir . "/${prefix}_DESeq2_volcanoEnhanced.png" );
    push( @result_files, $result_dir . "/${prefix}_geneAll_DESeq2-vsd-heatmap.png" );
    push( @result_files, $result_dir . "/${prefix}_geneAll_DESeq2-vsd-pca.png" );
    if ($exportSignificantGeneName) {
      push( @result_files, $result_dir . "/${prefix}_DESeq2_sig_genename.txt" );
    }
    if ($performWilcox) {
      push( @result_files, $result_dir . "/${prefix}_quantile_wilcox.csv" );
      push( @result_files, $result_dir . "/${prefix}_quantile_wilcox_sig.csv" );
    }
    push( @result_files, $result_dir . "/${comparison_name}.design" );

    #put volcanoPlot at end of result list for checking if the task has been finished.
    push( @result_files, $result_dir . "/${prefix}_DESeq2_volcanoPlot.png" );

    my $filtered = filter_array( \@result_files, $pattern, $removeEmpty );
    if ( scalar(@$filtered) > 0 || !$removeEmpty ) {
      $result->{$comparison_name} = $filtered;
    }
  }

  return $result;
}

1;
