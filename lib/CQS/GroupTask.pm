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

1;
