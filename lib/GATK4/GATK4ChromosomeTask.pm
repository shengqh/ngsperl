#!/usr/bin/perl
package GATK4::GATK4ChromosomeTask;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_docker_prefix} = "gatk4_";
  $self->{_export_home} = 1;
  bless $self, $class;
  return $self;
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $chromosomeStr = get_option($config, $section, "chromosome_names");
  my @chromosomes = split /,/, $chromosomeStr;
  
  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $result = {};

  for my $chr (@chromosomes) {
    my $chrTaskName = $task_name . "." . $chr;
    $result->{$chrTaskName} = $self->get_pbs_filename( $pbs_dir, $chrTaskName );
  }

  return($result);
}

1;
