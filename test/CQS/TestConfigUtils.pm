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
for my $key (keys %$unsorted){
  print $key;
}

my $sorted = get_raw_files($def, "sorted_section");
for my $key (keys %$sorted){
  print $key;
}

1;
