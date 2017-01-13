#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 2;
use CQS::StringUtils;

is_deeply( string_combination( [[ 'a', 'b', '' ], [ '1', '2', '' ]], ':' ), ["a:1", "a:2", "a", "b:1", "b:2", "b", "1", "2", "" ] );

is_deeply( string_repeat( [ 'a', 'b' ], 2 ), ["a", "a", "b", "b" ] );

1;
