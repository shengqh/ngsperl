#!/usr/bin/perl
package GATK4::IntervalsScatterTask;

use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use Tie::IxHash;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::ScatterTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::ScatterTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = __PACKAGE__;
  $self->{_docker_prefix} = "gatk4_";
  bless $self, $class;
  return $self;
}

#return a list of string
sub get_scatter_names {
  my ($self, $config, $section) = @_;
  my $interval_files = get_interval_file_map($config, $section);
  my $result = [keys %$interval_files];
  return ($result);
}

1;
