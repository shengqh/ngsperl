#!/usr/bin/perl
package CQS::TestGroupPickTask;

use strict;
use warnings;
use CQS::GroupPickTask;
use Test::More tests => 2;
use Data::Dumper;

my $test = CQS::GroupPickTask->new();

my $config = {
  "test" => {
    source => {
      "S1_1" => ["S1_1.txt"],
      "S1_2" => ["S1_2.txt"],
      "S2_1" => ["S2_1.txt"],
      "S2_2" => ["S2_2.txt"],
      "S3_1" => ["S3_1.txt"],
    },
    groups => {
      "G1" => ["S1_1", "S1_2"],
      "G2" => ["S2_1", "S2_2"],
      "G3" => ["S3_1"],
    },
    sample_index_in_group => 0,
  }
};

my $picked_0 = $test->result($config, "test");

#print(Dumper(%$source));

is_deeply($picked_0, 
  { "G1" => ["S1_1.txt"],
    "G2" => ["S2_1.txt"],
    "G3" => ["S3_1.txt"],
  });

$config->{"test"}{sample_index_in_group} = 1;
my $picked_1 = $test->result($config, "test");

#print(Dumper(%$source));

is_deeply($picked_1, 
  { "G1" => ["S1_2.txt"],
    "G2" => ["S2_2.txt"]
  });

1;


1
