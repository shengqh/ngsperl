#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 1;
use scRNA::Modules;

{#get_gmm_demux_option_map
  my $def = {
    HTO_samples => {
      "S1" => {
        "C1" => "PBMC_1975",
        "C2" => "PBMC_1932",
      },
      "S2" => {
        "T1" => "Aorta_1727",
        "T2" => "Aorta_1975",
      },
    }
  };
  my $map = get_gmm_demux_option_map($def);
  is_deeply($map, {
    "S1" => "C1,C2",
    "S2" => "T1,T2",
  });
}

1;