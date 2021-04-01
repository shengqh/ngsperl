#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use Test::More tests => 1;
use Pipeline::CRISPRScreen;

{
  my $batch_values= {
    'S1_1' => 1,
    'S1_2' => 1,
    'S2_1' => 2,
    'S2_2' => 2,
  };

  writeBatchFile($batch_values, '/scratch/cqs/shengq2/temp/test_batch.txt');
}


{
  my $batch_values= {
    'S1_1' => {
      Batch => 1,
      Condition => 1,
    },
    'S1_2' => {
      Batch => 1,
      Condition => 2,
    },
    'S2_1' => {
      Batch => 2,
      Condition => 1,
    },
    'S2_2' => {
      Batch => 2,
      Condition => 2
    },
  };

  writeBatchFile($batch_values, '/scratch/cqs/shengq2/temp/test_batch.txt');
}

1;
