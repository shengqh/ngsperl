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

#Return 
#{
#  groupName1 => [
#    [sampleName1_1, sampleFile1_1_1, sampleFile1_1_2],
#    [sampleName1_2, sampleFile1_2_1, sampleFile1_2_2],
#  ],
#  groupName2 => [
#    [sampleName2_1, sampleFile2_1_1, sampleFile2_1_2],
#    [sampleName2_2, sampleFile2_2_1, sampleFile2_2_2],
#  ],
#}
sub get_group_sample_map {
  my ( $self, $config, $section, $samplePattern ) = @_;

  my $rawFiles = get_raw_files( $config, $section, "source", $samplePattern );
  my $groups = get_raw_files( $config, $section, "groups" );
  my %group_sample_map = ();
  for my $groupName ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$groupName} };
    my @gfiles  = ();
    foreach my $sampleName (@samples) {
      my @bamFiles = @{ $rawFiles->{$sampleName} };
      my @sambam = ( $sampleName, @bamFiles );
      push( @gfiles, \@sambam );
    }
    $group_sample_map{$groupName} = \@gfiles;
  }

  return \%group_sample_map;
}

#Return 
#{
#  groupName1 => [sampleFile1_1_1, sampleFile1_1_2, sampleFile1_2_1, sampleFile1_2_2],
#  groupName2 => [sampleFile2_1_1, sampleFile2_1_2, sampleFile2_2_1, sampleFile2_2_2],
#}
sub get_group_samplefile_map {
  my ( $self, $config, $section, $samplePattern ) = @_;

  my $rawFiles = get_raw_files( $config, $section, "source", $samplePattern );
  my $groups = get_raw_files( $config, $section, "groups" );
  my %group_sample_map = ();
  for my $groupName ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$groupName} };
    my @gfiles  = ();
    foreach my $sampleName (@samples) {
      my @bamFiles = @{ $rawFiles->{$sampleName} };
      foreach my $bamFile (@bamFiles) {
        push( @gfiles, $bamFile );
      }
    }
    $group_sample_map{$groupName} = \@gfiles;
  }

  return \%group_sample_map;
}

1;
