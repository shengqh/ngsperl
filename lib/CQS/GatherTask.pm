#!/usr/bin/perl
package CQS::GatherTask;

use strict;
use warnings;
use File::Basename;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::Task;
use Data::Dumper;

our @ISA = qw(CQS::Task);

#This class is for task which generates one pbs for multiple input file.
#For example, combine count file to count table.
#For example, merge Impute2 result
sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_gt";
  $self->{_pbskey} = "", bless $self, $class;
  return $self;
}

#return a gather_name/scatter_name_list map
sub get_gather_map {
  my ( $self, $config, $section ) = @_;
  die "Override get_gather_map of " . $self->{_name} . " first.";
}

sub get_result_files {
  my ( $self, $config, $section, $result_dir, $gather_name ) = @_;

  die "Override get_result_files of " . $self->{_name} . " first.";
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $result = {};

  my $gather_maps = $self->get_gather_map( $config, $section );
  for my $gather_name ( keys %$gather_maps ) {
    $result->{$gather_name} = $self->get_pbs_filename( $pbs_dir, $gather_name );
  }

  return $result;
}

sub get_pbs_source {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $result = {};

  my $raw_files  = get_raw_files( $config, $section );
  my $gather_map = $self->get_gather_map( $config, $section );
  for my $gather_name ( keys %$gather_map ) {
    my $scatter_names = $gather_map->{$gather_name};
    if ( $self->acceptSample( $config, $section, $gather_name ) ) {
      my $pbs_file = $self->get_pbs_filename( $pbs_dir, $gather_name );
      $result->{$pbs_file} = $scatter_names;
    }
  }

  return $result;
}

sub result {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $result = {};

  my $raw_files  = get_raw_files( $config, $section );
  my $gather_map = $self->get_gather_map( $config, $section );
  for my $gather_name ( sort keys %$gather_map ) {
    if ( $self->acceptSample( $config, $section, $gather_name ) ) {
      $result->{$gather_name} = $self->get_result_files( $config, $section, $result_dir, $gather_name );
    }
  }

  return $result;
}

sub get_result_pbs {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $resultFiles = $self->result( $config, $section );
  my $pbsFiles    = $self->get_pbs_files( $config, $section );

  my $result = {};

  for my $gather_name ( sort keys %$resultFiles ) {
    $result->{$gather_name} = $pbsFiles->{$gather_name};
  }

  return $result;
}

1;
