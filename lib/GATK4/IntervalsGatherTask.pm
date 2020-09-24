#!/usr/bin/perl
package GATK4::IntervalsGatherTask;

use strict;
use warnings;
use File::Basename;
use List::Util qw[min];
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::GatherTask;
use GATK4::VariantFilterUtils;

our @ISA = qw(CQS::GatherTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_ig";
  $self->{_docker_prefix} = "gatk4_";
  $self->{_export_home} = 1;
  bless $self, $class;
  return $self;
}

#return a gather_name/scatter_name_list map
sub get_gather_map {
  my ($self, $config, $section) = @_;

  my $interval_files = get_interval_file_map($config, $section);
  my $raw_files = get_raw_files($config, $section);
  
  my $result = {};

  for my $sample_name (sort keys %$raw_files) {
    my $key_names = [];
    for my $scatter_name (sort keys %$interval_files) {
      my $key_name = get_key_name($sample_name, $scatter_name);
      push(@$key_names, $key_name);
    }
    $result->{$sample_name} = $key_names;
  }
  
  return $result;  
}

1;
