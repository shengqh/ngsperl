#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 2;
use Utils::FileUtils;

is( changeExtension( "a.b.c", ".d" ), "a.b.d", "changeExtension with extension" );
is( changeExtension( "a",     ".d" ), "a.d",   "changeExtension without extension" );

1;
