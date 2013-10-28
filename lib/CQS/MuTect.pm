#!/usr/bin/perl
package CQS::MuTect;

use strict;
use warnings;

sub new {
  require "GATK/MuTect.pm";
  my $class = "GATK::Mutect";
  return $class->new();
}

1;
