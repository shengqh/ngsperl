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
  bless $self, $class;
  return $self;
}

#return a gather_name/scatter_name_list map
sub get_gather_map {
  my ($self, $config, $section) = @_;

  my $interval_files = get_interval_file_map($config, $section);
  my $raw_files = get_raw_files($config, $section, "gather_name");
  
  my $result = {};

  for my $gather_name (sort keys %$raw_files) {
    my $key_names = [];
    for my $scatter_name (sort keys %$interval_files) {
      my $key_name = get_key_name($gather_name, $scatter_name);
      push(@$key_names, $key_name);
    }
    $result->{$gather_name} = $key_names;
  }
  
  return $result;  
}

sub get_gather_file_map {
  my ($self, $gather_map, $raw_files) = @_;

  my $result = {};
  for my $gather_name (sort keys %$gather_map) {
    my $key_names = $gather_map->{$gather_name};
    my $key_files = [];
    for my $key_name (@$key_names){
      my $cur_files = $raw_files->{$key_name};
      for my $cur_file (@$cur_files){
        push(@$key_files, $cur_file);
      }
    }
    $result->{$gather_name} = $key_files;
  }

  return($result);
}

1;
