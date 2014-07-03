#!/usr/bin/perl
package CQS::GroupTask;

use strict;
use warnings;
use CQS::Task;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_pbskey} = "groups";
  bless $self, $class;
  return $self;
}

sub get_group_sample_map {
  my ( $self, $config, $section ) = @_;

  my $rawFiles = get_raw_files( $config, $section );
  my $groups = get_raw_files( $config, $section, "groups" );
  my %group_sample_map = ();
  for my $groupName ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$groupName} };
    my @gfiles  = ();
    foreach my $sampleName (@samples) {
      my @bamFiles = @{ $rawFiles->{$sampleName} };
      my @sambam = ($sampleName, $bamFiles[0]);
      push( @gfiles, \@sambam );
    }
    $group_sample_map{$groupName} = \@gfiles;
  }

  return \%group_sample_map;
}

1;
