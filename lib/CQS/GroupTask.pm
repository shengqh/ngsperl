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
  $self->{_depend_all} = 1;
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

sub get_pbs_source {
  my ( $self, $config, $section ) = @_;

  my $pbsFiles = $self->get_pbs_files( $config, $section );
  my $result = {};
  my $groups = $config->{$section}{"groups"};
  if (defined $groups){
    for my $resKey ( keys %$pbsFiles ) {
      $result->{ $pbsFiles->{$resKey}} = $groups->{$resKey};
    }
  }else{
    for my $resKey ( keys %$pbsFiles ) {
      $result->{ $pbsFiles->{$resKey} } = [$resKey];
    }
  }
  return $result;
}

1;
