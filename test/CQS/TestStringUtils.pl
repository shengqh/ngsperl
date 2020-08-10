#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 3;
use CQS::StringUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

is_deeply( string_combination( [[ 'a', 'b', '' ], [ '1', '2', '' ]], ':' ), ["a:1", "a:2", "a", "b:1", "b:2", "b", "1", "2", "" ] );

is_deeply( string_repeat( [ 'a', 'b' ], 2 ), ["a", "a", "b", "b" ] );

my $old = {
  a => {
    a1 => {
      a11 => "old a11",
      a12 => "old a12",
    },
    a2 => "unchanged a2",
  },
  b => ["old b1","old b2","old b3"],
  c => "unchanged c"
};

my $new = {
  a => {
    a1 => {
      a11 => "new a11",
      a12 => "new a12",
    },
  },
  b => ["new b1", "new b2"],
};

my $res = merge($old, $new);

my $expect = {
          'b' => [
            "old b1","old b2","old b3",
                   'new b1',
                   'new b2'
                 ],
          'c' => 'unchanged c',
          'a' => {
                   'a1' => {
                             'a12' => 'new a12',
                             'a11' => 'new a11'
                           },
                   'a2' => 'unchanged a2'
                 }
        };
#print(Dumper($res));
is_deeply($res, $expect);

1;
