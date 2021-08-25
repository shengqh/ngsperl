#!/usr/bin/perl
package CQS::TestSamplePickTask;

use strict;
use warnings;
use CQS::SamplePickTask;
use Test::More tests => 2;
use Data::Dumper;

my $test = CQS::SamplePickTask->new();

my $config = {
  "test" => {
    source => {
      "S1_1" => ["S1_1.txt"],
      "S1_2" => ["S1_2.txt"],
      "S2_1" => ["S2_1.txt"],
      "S2_2" => ["S2_2.txt"],
      "S3_1" => ["S3_1.txt"],
    },
    sample_names => ["S1_1", "S1_2"],
  }
};

my $picked_0 = $test->result($config, "test");

#print(Dumper(%$source));

is_deeply($picked_0, 
  { 
    "S1_1" => ["S1_1.txt"],
    "S1_2" => ["S1_2.txt"],
  });

delete $config->{test}{sample_names};
$config->{test}{not_sample_names} = ["S1_1", "S1_2"];

my $picked_not = $test->result($config, "test");

#print(Dumper($picked_not));

is_deeply($picked_not, 
  { 
    "S2_1" => ["S2_1.txt"],
    "S2_2" => ["S2_2.txt"],
    "S3_1" => ["S3_1.txt"],
  });

1
