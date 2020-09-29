#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 11;
use CQS::ConfigUtils;
use Data::Dumper;

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

ok( option_contains_arg( "-i __SAMPLE__",              "-i" ) );
ok( option_contains_arg( "-o __OUTPUT__ -i",           "-i" ) );
ok( option_contains_arg( "-o __OUTPUT__ -i __INPUT__", "-i" ) );
ok( !option_contains_arg( "-o __OUTPUT__ -impossible __INPUT__", "-i" ) );

{    #test read_table
  my ( $tbl, $names ) = read_table( "../../data/cnv.txt", 3 );
  ok( 9 == scalar( keys %$tbl ) );
  is_deeply(
    $tbl->{'chr1:1705591-1705782'},
    {
      'Start'    => '1705342',
      'P_64B'    => '',
      'P_273_21' => '',
      'P_273_39' => '',
      'P_175_09' => '',
      'P_31B'    => '',
      'P_273_37' => 'DUP,2,4,33',
      'P_181'    => '',
      'P_56B'    => '',
      'P_196'    => '',
      'P_C04'    => 'DUP,2,4,30',
      'P_175_23' => 'DUP,2,4,44',
      'P_12B'    => '',
      'P_273_40' => '',
      'P_196F1'  => '',
      'P_175_18' => '',
      'P_175_10' => '',
      '#Chrom'   => 'chr1',
      'P_273_13' => '',
      'Gene'     => 'CDK11A',
      'P_42B'    => '',
      'P_273_22' => '',
      'P_196F2'  => '',
      'P_175_27' => 'DUP,2,4,32',
      'P_66B'    => '',
      'P_175_19' => 'DUP,2,4,42',
      'P_273_09' => 'DUP,2,4,48',
      'P_23B'    => '',
      'P_67B'    => '',
      'P_181F1'  => '',
      'P_175_12' => '',
      'P_273_38' => '',
      'P_273_15' => '',
      'P_175_33' => 'DUP,2,4,36',
      'P_175_06' => '',
      'End'      => '1706032',
      'P_38B'    => '',
      'P_C08'    => '',
      'P_273_20' => '',
      'P_60B'    => '',
      'P_09B'    => '',
      'P_273_03' => ''
    }
  );

  #print( Dumper($names) );
  is_deeply(
    $names,
    {
      '#Chrom'   => 1,
      'Start'    => 1,
      'End'      => 1,
      'Gene'     => 1,
      'P_09B'    => 1,
      'P_12B'    => 1,
      'P_175_06' => 1,
      'P_175_09' => 1,
      'P_175_10' => 1,
      'P_175_12' => 1,
      'P_175_18' => 1,
      'P_175_19' => 1,
      'P_175_23' => 1,
      'P_175_27' => 1,
      'P_175_33' => 1,
      'P_181'    => 1,
      'P_181F1'  => 1,
      'P_196'    => 1,
      'P_196F1'  => 1,
      'P_196F2'  => 1,
      'P_23B'    => 1,
      'P_273_03' => 1,
      'P_273_09' => 1,
      'P_273_13' => 1,
      'P_273_15' => 1,
      'P_273_20' => 1,
      'P_273_21' => 1,
      'P_273_22' => 1,
      'P_273_37' => 1,
      'P_273_38' => 1,
      'P_273_39' => 1,
      'P_273_40' => 1,
      'P_31B'    => 1,
      'P_38B'    => 1,
      'P_42B'    => 1,
      'P_56B'    => 1,
      'P_60B'    => 1,
      'P_64B'    => 1,
      'P_66B'    => 1,
      'P_67B'    => 1,
      'P_C04'    => 1,
      'P_C08'    => 1
    }
  );
}

1;
