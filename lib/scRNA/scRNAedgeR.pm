#!/usr/bin/perl
package scRNA::scRNAedgeR;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::UniqueR;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::UniqueR);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "";
  bless $self, $class;
  return $self;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $result = {};

  my $cluster_file = parse_param_file( $config, $section, "parameterFile1" );

  open( my $data, '<', $cluster_file ) or die "Could not open '$cluster_file' $!\n";
  my $line = <$data>;
  my %clusters;
  while ( $line = <$data> ) {
    chomp $line;
    my @fields = split ",", $line;
    $clusters{ $fields[3] } = 0;
  }
  my @uniqueClusters = sort keys %clusters;
  print "@uniqueClusters\n";

  my $pairs = get_raw_files( $config, $section, "parameterSampleFile2" );
  my @pairNames = sort keys %$pairs;
  print "@pairNames\n";

  for my $cluster (@uniqueClusters) {
    for my $pairName (@pairNames) {
      my $key    = $cluster . "." . $pairName;
      my $prefix = $result_dir . "/" . $task_name . "." . $key . ".edgeR";
      $result->{$key} = filter_array([ $prefix . ".csv", $prefix . "_sig_genename.txt", $prefix . "_GSEA.rnk" ], $pattern );
    }
  }

  return $result;
}

1;
