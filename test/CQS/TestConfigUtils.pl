#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 3;
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
    ".order" => [ "C", "B", "A" ],
    ".col"   => [ "CC", "BB", "AA" ]
  },
  unsorted_section => {
    source_ref => "unsorted"
  },
  sorted_section => {
    source_ref => "sorted"
  },
};

my $unsorted = get_raw_files( $def, "unsorted_section" );
my @unsortedKeys = keys %$unsorted;
is_deeply( \@unsortedKeys, [ "A", "B", "C" ] );

my $sorted = get_raw_files( $def, "sorted_section" );
my @sortedKeys = keys %$sorted;
is_deeply( \@sortedKeys, [ "C", "B", "A" ] );

my $sortedAttr = get_raw_files_attributes( $def, "sorted_section" );
is_deeply(
  $sortedAttr,
  {
    ".order" => [ "C",  "B",  "A" ],
    ".col"   => [ "CC", "BB", "AA" ]
  }
);

1;
