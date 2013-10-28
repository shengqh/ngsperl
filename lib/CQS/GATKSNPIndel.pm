#!/usr/bin/perl
package CQS::GATKSNPIndel;

use strict;
use warnings;

sub new {
  require "GATK/SNPIndel.pm";
  my $class = "GATK::SNPIndel";
  return $class->new();
}

1;
