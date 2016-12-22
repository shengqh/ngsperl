#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 2;
use CQS::ConfigUtils;

my $def = {
  unsorted => {
    B => 1,
    A => 1,
    C => 1
  },
  sorted => {
    B        => 1,
    A        => 1,
    C        => 1,
    ".order" => [ "C", "B", "A" ]
  },
  unsorted_section => {
    source_ref => "unsorted"
  },
  sorted_section => {
    source_ref => "sorted"
  },
};

my $unsorted = get_raw_files($def, "unsorted_section");
my @unsortedKeys = keys %$unsorted;
is_deeply(\@unsortedKeys, ["A", "B", "C"] );

my $sorted = get_raw_files($def, "sorted_section");
my @sortedKeys = keys %$sorted;
is_deeply(\@sortedKeys, ["C", "B", "A"] );

1;
