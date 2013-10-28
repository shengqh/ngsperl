#!/usr/bin/perl
package CQS::MuTect;

use strict;
use warnings;

sub new {
  require "GATK/Mutect.pm";
  my $class = "GATK::Mutect";
  return $class->new();
}

1;
