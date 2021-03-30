#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use Test::More tests => 1;
use Pipeline::CRISPRScreen;

{
  my $batch_values=[['S1_1', 'S1_2'], ['S2_1', 'S2_2']];

  writeBatchFile($batch_values, '/scratch/cqs/shengq2/temp/test_batch.txt');
}

1;
