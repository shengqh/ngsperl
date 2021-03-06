#!/usr/bin/perl
package GWAS::RiskScore;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::UniqueTask;
use CQS::NGSCommon;
use CQS::StringUtils;
use GWAS::GwasUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_pq";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $rawFiles = get_raw_files( $config, $section );
  my $scoreDefinitionFile = get_option_file( $config, $section, "score_definition_file" );

  for my $sampleName ( sort keys %$rawFiles ) {
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sampleName );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sampleName );
    my $log_desc = $cluster->get_log_description($log);

    my $files  = $rawFiles->{$sampleName};
    my $file   = $files->[0];
    my $prefix = getPlinkPrefix($file);

    my $final_file = $sampleName . ".profile";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
    print $pbs "
    plink --bfile $prefix --score $scoreDefinitionFile --out $sampleName
";
    $self->close_pbs( $pbs, $pbs_file );
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $rawFiles = get_raw_files( $config, $section );

  my $result = {};
  for my $sampleName ( sort keys %$rawFiles ) {
    $result->{$sampleName} = filter_array( ["$result_dir/${sampleName}.profile"], $pattern );
  }
  return $result;
}

1;
