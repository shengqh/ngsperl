#!/usr/bin/perl
package GATK4::GATK4Task;

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
  $self->{_name}   = __PACKAGE__;
  $self->{_docker_prefix} = "gatk4_";
  $self->{_export_home} = 1;
  bless $self, $class;
  return $self;
}

1;
