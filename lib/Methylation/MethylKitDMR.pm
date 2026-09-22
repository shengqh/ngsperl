#!/usr/bin/perl
package Methylation::MethylKitDMR;

use strict;
use warnings;
use CQS::ConfigUtils;
use CQS::IndividualR;

our @ISA = qw(CQS::IndividualR);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = __PACKAGE__;
  bless $self, $class;
  return $self;
}

sub get_pbs_source {
  my ( $self, $config, $section ) = @_;

  my $comparisons = get_raw_files( $config, $section, "source" );
  my $groups      = get_raw_files( $config, $section, "parameterSampleFile4" );
  my $pbs_files   = $self->get_pbs_files( $config, $section );
  my $result      = {};

  for my $comparison_name ( keys %$pbs_files ) {
    my ( $is_paired, $group_names ) = get_pair_groups( $comparisons, $comparison_name );
    die "Comparison $comparison_name must contain exactly two groups for methylKit DMR analysis.\n"
      if !defined($group_names) || scalar(@$group_names) != 2;

    my @samples = ();
    for my $group_name (@$group_names) {
      die "Cannot find group $group_name for methylKit DMR comparison $comparison_name.\n"
        if !defined $groups->{$group_name};
      die "Group $group_name has no samples for methylKit DMR comparison $comparison_name.\n"
        if scalar( @{ $groups->{$group_name} } ) == 0;
      push( @samples, @{ $groups->{$group_name} } );
    }
    $result->{ $pbs_files->{$comparison_name} } = \@samples;
  }

  return $result;
}

1;
