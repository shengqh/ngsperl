#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 3;
use CQS::StringUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

is_deeply( string_combination( [[ 'a', 'b', '' ], [ '1', '2', '' ]], ':' ), ["a:1", "a:2", "a", "b:1", "b:2", "b", "1", "2", "" ] );

is_deeply( string_repeat( [ 'a', 'b' ], 2 ), ["a", "a", "b", "b" ] );

my $res_p = capture_regex_groups("d_56kb_cMoP_H3K27ac_G2", "(^.+)(?:_input|_H3K27ac)(_.)");
ok($res_p eq "d_56kb_cMoP_G");

1;
