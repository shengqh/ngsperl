#!/usr/bin/perl
package CQS::FastqTrimmer;

use strict;
use warnings;
use Trimmer::RemoveFastqTerminalN;

our @ISA = qw(Trimmer::RemoveFastqTerminalN);

1;
