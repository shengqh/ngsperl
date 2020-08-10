#!/usr/bin/perl
package CQS::FilePickTask;

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
  my $pick_index = get_option( $config, $section, "sample_index", 0);

  my $result = {};

  for my $sample_name (keys %$raw_files) {
    my $sample_files = $raw_files->{$sample_name};
    my $result_files = [$sample_files->[$pick_index]];
    $result->{$sample_name} = filter_array( $result_files, $pattern );
  }

  return $result;
}

1;
