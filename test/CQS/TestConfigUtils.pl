#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 8;
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

my $def2 = {
    groups => {
      "g1" => [ "g1s1", "g1s2" ],
      "g2" => [ "g2s1", "g2s2" ],
      "g3" => [ "g3s1", "g3s2" ],
    },
    correlation_groups => {
      "g1g2" => {
        groups => [ "g1", "g2" ],
        cov    => [ 1,    2, 1, 2 ],
      },
      "g1g3" => [ "g1", "g3" ],
    },
};

my $map = get_pair_group_sample_map( $def2->{correlation_groups}, $def2->{groups} );
is_deeply(
  $map,
  {
    "g1g2" => {
      "g1" => [ "g1s1", "g1s2" ],
      "g2" => [ "g2s1", "g2s2" ],
    },
    "g1g3" => {
      "g1" => [ "g1s1", "g1s2" ],
      "g3" => [ "g3s1", "g3s2" ],
    },
  }
);

ok(option_contains_arg("-i __SAMPLE__", "-i"));
ok(option_contains_arg("-o __OUTPUT__ -i", "-i"));
ok(option_contains_arg("-o __OUTPUT__ -i __INPUT__", "-i"));
ok(! option_contains_arg("-o __OUTPUT__ -impossible __INPUT__", "-i"));

1;
