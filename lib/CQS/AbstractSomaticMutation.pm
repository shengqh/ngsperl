#!/usr/bin/perl
package CQS::AbstractSomaticMutation;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "AbstractSomaticMutation";
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
    my $index   = 0;
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
