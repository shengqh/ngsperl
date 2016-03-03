#!/usr/bin/perl
package CQS::GroupTask;

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
  $self->{_pbskey} = "groups";
  bless $self, $class;
  return $self;
}

sub get_grouped_raw_files {
  my ( $self, $config, $section, $group_key ) = @_;
  my $raw_files;
  if ( has_raw_files( $config, $section, $group_key ) ) {
    $raw_files = get_group_samplefile_map_key( $config, $section, "", $group_key );
  }
  else {
    $raw_files = get_raw_files( $config, $section );
  }
  return $raw_files;
}

1;
