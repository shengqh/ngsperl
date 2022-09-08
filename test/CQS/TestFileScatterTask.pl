#!/usr/bin/perl
package CQS::TestFileScatterTask;

use strict;
use warnings;
use CQS::FileScatterTask;
use Test::More tests => 3;
use Data::Dumper;

my $test = CQS::FileScatterTask->new();

my $config = {
  "test" => {
    source => {
      "S" => ["S1_1.txt", "S1_2.txt", "S2_1.txt", "S2_2.txt", "S3_1.txt", "S3_2.txt"],
    },
  }
};

my $picked = $test->result($config, "test");

#print(Dumper($picked));

is_deeply($picked, 
  { 
          'S_ITER_01' => [
                           'S1_1.txt'
                         ],
          'S_ITER_02' => [
                           'S1_2.txt'
                         ],
          'S_ITER_03' => [
                           'S2_1.txt'
                         ],
          'S_ITER_04' => [
                           'S2_2.txt'
                         ],
          'S_ITER_05' => [
                           'S3_1.txt'
                         ],
          'S_ITER_06' => [
                           'S3_2.txt'
                         ]
  });

$config->{"test"}{step} = 2;
$picked = $test->result($config, "test");

#print(Dumper($picked));

is_deeply($picked, 
  { 
          'S_ITER_01' => [
                           'S1_1.txt',
                           'S1_2.txt'
                         ],
          'S_ITER_02' => [
                           'S2_1.txt',
                           'S2_2.txt'
                         ],
          'S_ITER_03' => [
                           'S3_1.txt',
                           'S3_2.txt'
                         ]
  });

$config->{"test"}{step} = 3;
$picked = $test->result($config, "test");

#print(Dumper($picked));

is_deeply($picked, 
  { 
         'S_ITER_01' => [
                           'S1_1.txt',
                           'S1_2.txt',
                           'S2_1.txt'
                         ],
          'S_ITER_02' => [
                           'S2_2.txt',
                           'S3_1.txt',
                           'S3_2.txt'
                         ]
  });

1;
