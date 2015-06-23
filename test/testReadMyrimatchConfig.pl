#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Test::Simple;
use Proteomics::Engine::Myrimatch;

my $mm = Proteomics::Engine::Myrimatch->new();

my $hash = $mm->readConfigurationFile(dirname(__FILE__) . "/../data/myrimatch.cfg");

print "$_ => $hash->{$_}\n" for keys %{$hash};

1;
