#!/usr/bin/perl
package CQS::ScatterTask;

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

#This class is for task which generates multiple pbs/result files for each input file.
#For example, CombineGVCFs by chromosome in GATK4 SNV pipeline.
#For example, Impute2 which generate multiple pbs files and multiple result files for each sample
sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_sc";
  $self->{_pbskey} = "",
  bless $self, $class;
  return $self;
}

#get sample names, default is from "source"
sub get_sample_names {
  my ($self, $config, $section) = @_;
  my $raw_files = get_raw_files($config, $section);
  return [sort keys %$raw_files];
}

#return a list of string
sub get_scatter_names {
  my ($self, $config, $section) = @_;
  die "Override get_scatter_names of " . $self->{_name} . " first.";
}

#return result based on $result_dir, $sample_name, $key_name.
sub get_result_files {
  my ( $self, $config, $section, $result_dir, $sample_name, $scatter_name, $key_name ) = @_;

  die "Override get_result_files of " . $self->{_name} . " first.";
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $result = {};

  my $sample_names = $self->get_sample_names($config, $section);
  my $scatter_names = $self->get_scatter_names($config, $section);
  for my $sample_name (@$sample_names){
    if ( $self->acceptSample( $config, $section, $sample_name ) ) {
      for my $scatter_name ( @$scatter_names ) {
        my $key_name = get_key_name($sample_name, $scatter_name);
        $result->{$key_name} = $self->get_pbs_filename( $pbs_dir, $key_name );
      }
    }
  }

  return $result;
}

sub get_pbs_source {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $result = {};

  my $sample_names = $self->get_sample_names($config, $section);
  my $scatter_names = $self->get_scatter_names($config, $section);
  for my $sample_name (@$sample_names){
    if ( $self->acceptSample( $config, $section, $sample_name ) ) {
      for my $scatter_name ( @$scatter_names ) {
        my $key_name = get_key_name($sample_name, $scatter_name);
        my $pbs_file = $self->get_pbs_filename( $pbs_dir, $key_name );
        $result->{$pbs_file} = [$sample_name];
      }
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

  my $sample_names = $self->get_sample_names($config, $section);
  my $scatter_names = $self->get_scatter_names($config, $section);
  for my $sample_name (@$sample_names){
    if ( $self->acceptSample( $config, $section, $sample_name ) ) {
      for my $scatter_name ( @$scatter_names ) {
        my $key_name = get_key_name($sample_name, $scatter_name);
        $result->{$key_name} = $self->get_result_files( $config, $section, $result_dir, $sample_name, $scatter_name, $key_name );
      }
    }
  }

  return $result;
}

1;