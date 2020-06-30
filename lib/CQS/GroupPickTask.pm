#!/usr/bin/perl
package CQS::GroupPickTask;

use strict;
use warnings;
use CQS::Task;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  bless $self, $class;
  return $self;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my $raw_files = get_raw_files($config, $section);
  my $groups = get_raw_files($config, $section, "groups");
  my $pick_index = get_option( $config, $section, "sample_index_in_group", 0);

  my $result = {};

  for my $group_name (keys %$groups) {
    my $samples = $groups->{$group_name};
    my $sample_name = $samples->[$pick_index];
    my $result_files = $raw_files->{$sample_name};
    $result->{$group_name} = filter_array( $result_files, $pattern );
  }

  return $result;
}

1;
