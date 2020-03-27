#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::ProgramWrapperOneToOne;
use Test::More tests => 4;

my $test = CQS::ProgramWrapperOneToOne->new();

my $config = {
  general => {
    task_name => "one2one",
  },
  "test" => {
    class => "CQS::ProgramWrapperManyToOne",
    target_dir => "/test",
    source     => {
      "G1" => [ "S1_1", "S1_2" ],
      "G2" => [ "S2_1", "S2_2" ],
    },
    output_file_ext => ".1.csv",
    output_other_ext => ".2.csv",
    output_to_same_folder => 1,
    pbs => {
      email => "a\@b.c",
    }
  },
};

#test result
my $expect_result = {
  'G1' => ['/test/result/G1.1.csv','/test/result/G1.2.csv'],
  'G2' => ['/test/result/G2.1.csv','/test/result/G2.2.csv'],
};
my $actual_result = $test->result( $config, "test" );
is_deeply( $actual_result, $expect_result );

#test pbs
my $expect_pbs = {
  'G1' => '/test/pbs/G1_o2o.pbs',
  'G2' => '/test/pbs/G2_o2o.pbs'
};
my $actual_pbs = $test->get_pbs_files( $config, "test" );
is_deeply( $actual_pbs, $expect_pbs );

#test result_pbs
my $expect_pbs_sample_map = {
  '/test/pbs/G1_o2o.pbs' => ['G1'],
  '/test/pbs/G2_o2o.pbs' => ['G2']
};

my $pbs_sample_map = $test->get_pbs_source($config, "test");
is_deeply( $pbs_sample_map, $expect_pbs_sample_map );

#test result_pbs
my $expect_result_pbs_map = {
  'G1' => '/test/pbs/G1_o2o.pbs',
  'G2' => '/test/pbs/G2_o2o.pbs',
};

my $result_pbs_map = $test->get_result_pbs($config, "test");
is_deeply( $result_pbs_map, $expect_result_pbs_map );

1