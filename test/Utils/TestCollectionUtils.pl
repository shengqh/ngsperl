#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 3;
use Utils::CollectionUtils;

my $map = readDictionaryByIndex("noheader.dic");

is( $map->{"root"}, "Others" );
is( $map->{"Promicromonospora sp. 10.25-Bb"}, "Bacteria" );

my $map2 = readDictionaryByColumnName("noheader.dic", "root", "Others");
is( $map2->{"Promicromonospora sp. 10.25-Bb"}, "Bacteria" );

1;
