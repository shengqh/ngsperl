#!/usr/bin/perl
package CQS::TestOneToOneTask;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use MyTest::OneToOneTask;
use Test::More tests => 8;

my $test = MyTest::OneToOneTask->new();

my $config = {
  general => {
    task_name => "one2one",
  },
  "test" => {
    class => "MyTest::OneToOneTask",
    target_dir => "/test",
    source     => {
      "G1" => [ "S1_1", "S1_2" ],
      "G2" => [ "S2_1", "S2_2" ],
    },
    pbs => {
      email => "a\@b.c",
    }
  },
  "test2" => {
    class => "MyTest::OneToOneTask",
    target_dir => "/test2",
    source_ref => ["test"],
    pbs => {
      email => "a\@b.c",
    }
  },
};

#test result
my $expect_result = {
  'G1' => ['/test/result/G1.csv'],
  'G2' => ['/test/result/G2.csv']
};
my $actual_result = $test->result( $config, "test" );
is_deeply( $actual_result, $expect_result );

#test pbs
my $expect_pbs = {
  'G1' => '/test/pbs/G1_11.pbs',
  'G2' => '/test/pbs/G2_11.pbs'
};
my $actual_pbs = $test->get_pbs_files( $config, "test" );
is_deeply( $actual_pbs, $expect_pbs );

is($test->get_final_file($config, "test", "/temp"), "/test/result/G2.csv");
is($test->get_final_file($config, "test", "/test"), "result/G2.csv");
is($test->get_final_file($config, "test", "/test/"), "result/G2.csv");
is($test->get_final_file($config, "test", "/test", "G1"), "result/G1.csv");
is($test->get_final_file($config, "test", "/test/result", "G1"), "G1.csv");

my $expect2_pbs = {
  'G1' => '/test2/pbs/G1_11.pbs',
  'G2' => '/test2/pbs/G2_11.pbs'
};
my $actual2_pbs = $test->get_pbs_files( $config, "test2" );
is_deeply( $actual2_pbs, $expect2_pbs );


1