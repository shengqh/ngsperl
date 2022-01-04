#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use File::Basename;
use CQS::ConfigUtils;
use Data::Dumper;
use Test::More tests => 42;

{ #test is_string
  ok(is_string("string"));
  ok(! is_string([]));
  ok(! is_string({}));
}

{ #test is_array
  ok(! is_array("string"));
  ok(is_array([]));
  ok(! is_array({}));
}

{ #test is_not_array
  ok(is_not_array("string"));
  ok(!is_not_array([]));
  ok(is_not_array({}));
}

{ #test is_hash
  ok(! is_hash("string"));
  ok(! is_hash([]));
  ok(is_hash({}));
}

{ #test is_not_hash
  ok(is_not_hash("string"));
  ok(is_not_hash([]));
  ok(!is_not_hash({}));
}

{
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
}

{ # test_get_pair_group_sample_map 
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
}

{ # test option_contains_arg
  ok( option_contains_arg( "-i __SAMPLE__",              "-i" ) );
  ok( option_contains_arg( "-o __OUTPUT__ -i",           "-i" ) );
  ok( option_contains_arg( "-o __OUTPUT__ -i __INPUT__", "-i" ) );
  ok( !option_contains_arg( "-o __OUTPUT__ -impossible __INPUT__", "-i" ) );
}

{    #test read_table
  my $datafile = File::Spec->rel2abs(dirname(__FILE__) . "/../../data/cnv.txt");
  my ( $tbl, $names ) = read_table( $datafile , 3 );
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

my $def3 = {
  target_dir => "/scratch/cqs/shengq2/temp",
  files => {
    "MB02v1" => ["/data/stein_lab/mjo_sRNA_data/rawdata_4829/4829-JS-1_1_S01_L005_R1_001.fastq.gz"],
    "MB02v2" => ["/data/stein_lab/mjo_sRNA_data/rawdata_4829/4829-JS-1_2_S01_L005_R1_001.fastq.gz"],
    "posJFS" => ["/data/stein_lab/mjo_sRNA_data/rawdata_4829/4829-JS-1_35_S01_L005_R1_001.fastq.gz"],
    "negctrl" => ["/data/stein_lab/mjo_sRNA_data/rawdata_4829/4829-JS-1_36_S01_L005_R1_001.fastq.gz"],
  },

  covariance_patterns => {
    subject => {
      pattern => "(.*)v"
    }
  },
};

$def3->{groups_pattern} = "(v\\d)";
my $groups = get_groups_by_pattern($def3);
my $expect_groups = {
          'negctrl' => [
                         'negctrl'
                       ],
          'v2' => [
                    'MB02v2'
                  ],
          'posJFS' => [
                        'posJFS'
                      ],
          'v1' => [
                    'MB02v1'
                  ]
        };
is_deeply( $groups, $expect_groups );     

$def3->{groups_pattern} = {
  "v1" => "v1",
  "v2" => "v2",
  "other" => "pos|neg"
};
$groups = get_groups_by_pattern($def3);
$expect_groups = {
         'v1' => [
                    'MB02v1'
                  ],
          'v2' => [
                    'MB02v2'
                  ],
          'other' => [
                       'negctrl',
                       'posJFS'
                     ]
                  };
is_deeply( $groups, $expect_groups );     

my ($cov_map, $covariances, $samplenames) = get_covariances_by_pattern($def3);
my $cov_expect = {
          'subject' => {
                         'negctrl' => 'negctrl',
                         'posJFS' => 'posJFS',
                         'MB02v1' => 'MB02',
                         'MB02v2' => 'MB02'
                       }
        };
is_deeply( $cov_map, $cov_expect );     

{ # test get_output_ext_list
  my $config = {
    "test" => {
      output_file_ext      => ".final.rds;.aaa.csv",
      output_other_ext  => ".cluster.csv; .allmarkers.csv,.top10markers.csv ;_ur.html; ",
    }
  };
  my $exts = get_output_ext_list($config, "test");
  #print(Dumper($exts));
  is_deeply( $exts, 
        [
          '.final.rds',
          '.aaa.csv',
          '.cluster.csv',
          '.allmarkers.csv',
          '.top10markers.csv',
          '_ur.html'
        ] );    

  is(get_output_ext($config, "test"), ".final.rds");
}

{ #test get_hash_level2
  my $curated_gene_dotplot = {
    "Microphage" => {
      clusters => ["1","6","20"],
      genes => {
        IL_17A_pathway_1 => [qw(IL17A IL17RC)],
        IL_17A_pathway_2 => [qw(CXCL3 CSF3)],
      }
    }
  };

  my $clusters = get_hash_level2($curated_gene_dotplot, "clusters");
  is_deeply( $clusters, 
    {
      "Microphage" => ["1","6","20"]
    });
  my $genes = get_hash_level2($curated_gene_dotplot, "genes");
  is_deeply( $genes, 
    {
      "Microphage" => {
        IL_17A_pathway_1 => [qw(IL17A IL17RC)],
        IL_17A_pathway_2 => [qw(CXCL3 CSF3)],
      }
    });
}

{ # get_expanded_genes, input array
  my $gene_array = [qw(G1 G2)];
  my $actual = get_expanded_genes($gene_array);
  is_deeply( $actual, 
    {
      "interesting_genes" => {
        clusters => [],
        genes => $gene_array
      }
    });
}

{ # get_expanded_genes, input array and hash
  my $curated_gene_dotplot = {
    "GeneArray" => [qw(G1 G2)],
    "Microphage" => {
      clusters => ["1","6","20"],
      genes => {
        IL_17A_pathway_1 => [qw(IL17A IL17RC)],
        IL_17A_pathway_2 => [qw(CXCL3 CSF3)],
      }
    }
  };

  my $actual = get_expanded_genes($curated_gene_dotplot);
  is_deeply( $actual, 
    {
      "GeneArray" => {
        clusters => [],
        genes =>  [qw(G1 G2)]
      },
      "Microphage" => {
        clusters => ["1","6","20"],
        genes => {
          IL_17A_pathway_1 => [qw(IL17A IL17RC)],
          IL_17A_pathway_2 => [qw(CXCL3 CSF3)],
        }
      }      
    });
}

{ # test parse_curated_genes
  my $curated_gene_dotplot = {
    "GeneArray" => [qw(G1 G2)],
    "Microphage" => {
      clusters => ["1","6","20"],
      genes => {
        IL_17A_pathway_1 => [qw(IL17A IL17RC)],
        IL_17A_pathway_2 => [qw(CXCL3 CSF3)],
      }
    }
  };

  my ($parsed, $clusters, $genes) = parse_curated_genes($curated_gene_dotplot);

  is_deeply( $clusters, 
    {
      "GeneArray" => [],
      "Microphage" => ["1","6","20"]
    });

  is_deeply( $genes, 
    {
      "GeneArray" => [qw(G1 G2)],
      "Microphage" => {
        IL_17A_pathway_1 => [qw(IL17A IL17RC)],
        IL_17A_pathway_2 => [qw(CXCL3 CSF3)],
      }
    });
}

{#test get_first_result_file
  my $config = {
    my_file1 => "file1",
    my_file2 => ["file2", "file3"],
    my_file3 => {task_name => ["file3", "file4"]},
  };

  is(get_first_result_file($config, "my_file1", ""),
    "file1");

  is(get_first_result_file($config, "my_file2", ""),
    "file2");

  is(get_first_result_file($config, "my_file3", ""),
    "file3");
}

{# test get_correlation_groups_by_pattern
  my $def = {
    files => {                                       
      'LDLR_PEL_01' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_1_S1_L005_R1_001.fastq.gz' ],
      'LDLR_PEL_02' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_2_S1_L005_R1_001.fastq.gz' ],
      'LDLR_PEL_03' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_3_S1_L005_R1_001.fastq.gz' ],
      'DKO_PEL_04' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_4_S1_L005_R1_001.fastq.gz' ],
      'DKO_PEL_05' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_5_S1_L005_R1_001.fastq.gz' ],
      'DKO_PEL_06' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_6_S1_L005_R1_001.fastq.gz' ],
      'LDLR_PEL_07' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_7_S1_L005_R1_001.fastq.gz' ],
      'LDLR_PEL_08' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_8_S1_L005_R1_001.fastq.gz' ],
      'LDLR_PEL_09' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_9_S1_L005_R1_001.fastq.gz' ],
      'DKO_PEL_10' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_10_S1_L005_R1_001.fastq.gz' ],
      'DKO_PEL_11' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_11_S1_L005_R1_001.fastq.gz' ],
      'DKO_PEL_12' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_12_S1_L005_R1_001.fastq.gz' ],
      'LDLR_SUP_13' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_13_S1_L005_R1_001.fastq.gz' ],
      'LDLR_SUP_14' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_14_S1_L005_R1_001.fastq.gz' ],
      'LDLR_SUP_15' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_15_S1_L005_R1_001.fastq.gz' ],
      'DKO_SUP_16' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_16_S1_L005_R1_001.fastq.gz' ],
      'DKO_SUP_17' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_17_S1_L005_R1_001.fastq.gz' ],
      'DKO_SUP_18' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_18_S1_L005_R1_001.fastq.gz' ],
      'LDLR_SUP_19' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_19_S1_L005_R1_001.fastq.gz' ],
      'LDLR_SUP_20' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_20_S1_L005_R1_001.fastq.gz' ],
      'LDLR_SUP_21' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_21_S1_L005_R1_001.fastq.gz' ],
      'DKO_SUP_22' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_22_S1_L005_R1_001.fastq.gz' ],
      'DKO_SUP_23' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_23_S1_L005_R1_001.fastq.gz' ],
      'DKO_SUP_24' => [ '/data/vickers_lab/20211227_7377_RA_NextFlex_smRNA_mouse/7377-RA-1_24_S1_L005_R1_001.fastq.gz' ],
    },
    
    groups_pattern => '(.+)_\d',
    
    "pairs" => {
      "PEL_DKO_vs_LDLR" => ["LDLR_PEL", "DKO_PEL"],
      "SUP_DKO_vs_LDLR" => ["LDLR_SUP", "DKO_SUP"],
      "LDLR_PEL_vs_SUP" => ["LDLR_SUP", "LDLR_PEL"],
      "DKO_PEL_vs_SUP" => ["DKO_SUP", "DKO_PEL"],
    },
  };

  $def->{groups} = get_groups_by_pattern($def);
  is_deeply($def->{groups} ,{
          'DKO_PEL' => [
                         'DKO_PEL_04',
                         'DKO_PEL_05',
                         'DKO_PEL_06',
                         'DKO_PEL_10',
                         'DKO_PEL_11',
                         'DKO_PEL_12'
                       ],
          'LDLR_SUP' => [
                          'LDLR_SUP_13',
                          'LDLR_SUP_14',
                          'LDLR_SUP_15',
                          'LDLR_SUP_19',
                          'LDLR_SUP_20',
                          'LDLR_SUP_21'
                        ],
          'DKO_SUP' => [
                         'DKO_SUP_16',
                         'DKO_SUP_17',
                         'DKO_SUP_18',
                         'DKO_SUP_22',
                         'DKO_SUP_23',
                         'DKO_SUP_24'
                       ],
          'LDLR_PEL' => [
                          'LDLR_PEL_01',
                          'LDLR_PEL_02',
                          'LDLR_PEL_03',
                          'LDLR_PEL_07',
                          'LDLR_PEL_08',
                          'LDLR_PEL_09'
                        ]
        });

  my $corr_groups = get_correlation_groups_by_pattern($def);

  is_deeply( $corr_groups, 
    {
      "all" => ["DKO_PEL", "DKO_SUP", "LDLR_PEL", "LDLR_SUP"],
    });
}

1;
