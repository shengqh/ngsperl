#!/usr/bin/perl
use strict;
use warnings;

use CQS::PBS;
use Test::Simple;

my $pbs8 = {
  "email"    => "email",
  "nodes"    => "1:ppn=8",
  "walltime" => "240",
  "mem"      => "40gb"
};

my $pbs1 = {
  "email"    => "email",
  "nodes"    => "1:ppn=1",
  "walltime" => "240",
  "mem"      => "40gb"
};

my ( $thread8 ) = get_pbs_thread( $pbs8 );
ok( $thread8 eq 8 );

my ( $thread1 ) = get_pbs_thread( $pbs1 );
ok( $thread1 eq 1 );

1;
