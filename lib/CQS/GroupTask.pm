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

sub get_pbs_key {
  my ($self, $config, $section) = @_;
  if ( has_raw_files( $config, $section, "groups" ) ) {
    return "groups";
  }else{
    return "source";
  }
}

1;
