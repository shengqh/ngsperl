#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 2;
use CQS::ConfigUtils;

my $config = {
  "test" => {
    "option-a" => 0.1,
    "option-c" => "aaa",
  }
};

is( get_parameter_options( $config, "test", "--", ["option-a", "option-b"] ), " --option-a 0.1", "get_parameter_options of ConfigUtils" );
is( get_parameter_options( $config, "test", "--", ["option-a", "option-c"] ), " --option-a 0.1 --option-c aaa", "get_parameter_options of ConfigUtils" );

1;
