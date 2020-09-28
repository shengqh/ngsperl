#!/usr/bin/perl
package GATK4::GATK4ChromosomeTask;

use strict;
use warnings;
use File::Basename;
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
  my $chromosomeStr = get_option($config, $section, "chromosome_names");
  my @chromosomes = split /,/, $chromosomeStr;
  return \@chromosomes;
}

1;
