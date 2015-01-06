#!/usr/bin/perl
package CQS::ClusterScript;

use strict;
use warnings;
use CQS::ConfigUtils;

sub new {
  my ($class) = @_;
  my $self = { _name => undef };
  bless $self, $class;
  return $self;
}

sub name {
  my ($self) = @_;
  return $self->{_name};
}

sub get_cluster_desc {
}

sub get_log_desc {
}

sub get_submit_command {
}

1;
