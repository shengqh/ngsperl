#!/usr/bin/perl
package Chipseq::AbstractMACSCallpeak;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::GroupTask;
use CQS::NGSCommon;
use CQS::StringUtils;
use Data::Dumper;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  bless $self, $class;
  return $self;
}

sub get_current_raw_files {
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
