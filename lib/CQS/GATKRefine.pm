#!/usr/bin/perl
package CQS::GATKRefine;

use strict;
use warnings;

sub new {
  require "GATK/Refine.pm";
  my $class = "GATK::Refine";
  return $class->new();
}

1;
