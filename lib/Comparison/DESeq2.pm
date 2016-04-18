#!/usr/bin/perl
package Comparison::DESeq2;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
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

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $comparisons = get_raw_files( $config, $section );
  my @comparison_names = sort keys %{$comparisons};

  my $totalPair = scalar(@comparison_names);
  if ( 0 == $totalPair ) {
    die "No pair defined!";
  }

  my $groups = get_raw_files( $config, $section, "groups" );

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

  my %tpgroups = ();
  for my $group_name ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$group_name} };
    $tpgroups{$group_name} = "\"" . join( "\",\"", @samples ) . "\"";
  }

  my $rfile = $result_dir . "/${task_name}.r";
  open( my $rf, ">$rfile" )     or die "Cannot create $rfile";
  open( my $rt, "<$rtemplate" ) or die $!;

  my $readfunc;
  my $readparam;
  if ( $countfile =~ /csv$/ ) {
    $readfunc  = "read.csv";
    $readparam = "";
  }
  else {
    $readfunc  = "read.table";
    $readparam = ", sep=\"\\t\"";
  }
  print $rf "
setwd(\"$result_dir\")  
  
data<-${readfunc}(\"$countfile\",row.names=1, header=T, check.names=F $readparam)

taskName<-'$task_name'
showLabelInPCA<-$showLabelInPCA
showDEGeneCluster<-$showDEGeneCluster
pvalue<-$pvalue
foldChange<-$foldChange
minMedianInGroup<-$minMedianInGroup

comparisons=list(";
  my $first = 0;
  for my $comparison_name (@comparison_names) {
    $first++;

    my $covariances = {};

    my $gNames = $comparisons->{$comparison_name};
    my @group_names;

    if ( ref $gNames eq ref {} ) {
      @group_names = @{ $gNames->{groups} };
      for my $key ( sort keys %$gNames ) {
        next if ( $key == "groups" );
        $covariances->{$key} = $gNames->{$key};
      }
    }
    else {
      @group_names = @{$gNames};
    }
    my @covariances_keys = sort keys %$covariances;

    #print( Dumper(@group_names) );

    if ( scalar(@group_names) != 2 ) {
      die "Comparison of $comparison_name should contains and only contains two groups!";
    }

    my $g1 = $group_names[0];
    my $g2 = $group_names[1];
    my @s1 = @{ $groups->{$g1} };
    my @s2 = @{ $groups->{$g2} };

    my $total_sample_count = scalar(@s1) + scalar(@s2);

    for my $key ( keys %$covariances ) {
      my $values = $covariances->{$key};
      if ( !( ref $values eq ref [] ) ) {
        die "Covariances of " . $key . " shoud be array reference!";
      }

      if ( scalar(@$values) != $total_sample_count ) {
        die "Number of covariance value of " . $key . " shoud be $total_sample_count !";
      }
    }

    my $filename = "${comparison_name}.design";
    if ( $first != 1 ) {
      print $rf ",";
    }
    print $rf "
  \"${comparison_name}\" = c(\"$filename\", \"$g1\", \"$g2\")";

    my $cdfile = $result_dir . "/$filename";
    open( my $cd, ">$cdfile" ) or die "Cannot create $cdfile";
    if ( scalar(@covariances_keys) > 0 ) {
      print $cd "Sample\tCondition\t", join( "\t", @covariances_keys ), "\n";
    }
    else {
      print $cd "Sample\tCondition\n";
    }
    for my $i ( 0 .. $#s1 ) {
      my $sname = $s1[$i];
      print $cd "${sname}\t${g1}";
      if ( scalar(@covariances_keys) > 0 ) {
        for my $key (@covariances_keys) {
          print $cd "\t" . $covariances->{$key}[$i];
        }
      }
      print $cd "\n";
    }
    for my $i ( 0 .. $#s2 ) {
      my $sname = $s2[$i];
      print $cd "${sname}\t${g2}";
      if ( scalar(@covariances_keys) > 0 ) {
        for my $key (@covariances_keys) {
          print $cd "\t" . $covariances->{$key}[ $i + scalar(@s1) ];
        }
      }
      print $cd "\n";
    }
    close $cd;
  }
  print $rf "
) \n\n";

  while (<$rt>) {
    if ( $_ =~ '^#' ) {
      next;
    }
    last;
  }
  while (<$rt>) {
    chomp;
    print $rf $_, "\n";
  }
  close $rt;
  close $rf;

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  my $final_file = $comparison_names[0] . "_min${minMedianInGroup}_DESeq2_sig.csv";
  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  print $pbs "R --vanilla -f $rfile \n";

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $comparisons = get_raw_files( $config, $section );
  my $minMedianInGroup = get_option( $config, $section, "min_median_read", 0 );

  my $result = {};
  for my $comparison_name ( sort keys %{$comparisons} ) {
    my @result_files = ();
    push( @result_files, $result_dir . "/${comparison_name}_min${minMedianInGroup}.csv" );
    push( @result_files, $result_dir . "/${comparison_name}_min${minMedianInGroup}_DESeq2.csv" );
    push( @result_files, $result_dir . "/${comparison_name}_min${minMedianInGroup}_DESeq2_sig.csv" );
    $result->{$comparison_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
