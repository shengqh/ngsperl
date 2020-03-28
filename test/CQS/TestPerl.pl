#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::Perl;
use Test::More tests => 4;

my $config = {
  general => {
    task_name => "test",
    cluster   => "slurm",
  },
  "task" => {
    class      => "CQS::Perl",
    perform    => 1,
    target_dir => "target",
    option     => "",
    perlFile   => "../SmallRNA/filterTrnaXml.py",
    source     => {
      "test1" => [ "test1_1.zip", "test1_2.zip" ],
      "test2" => [ "test2_1.zip", "test2_2.zip" ],
    },
    output_to_same_folder => 1,
    output_ext            => ".png",
    output_other_ext      => ".csv",
    sh_direct             => 1,
    pbs                   => {
      "email"     => "quanhu.sheng.1\@vumc.org",
      "emailType" => "FAIL",
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  },
};

my $test = CQS::Perl->new();

#test result
my $expect_result = {
  'test1' => [ 'target/result/test1.png', 'target/result/test1.csv' ],
  'test2' => [ 'target/result/test2.png', 'target/result/test2.csv' ],
};
my $actual_result = $test->result( $config, "task" );
is_deeply( $actual_result, $expect_result );

#test pbs
my $expect_pbs = {
  'test1' => 'target/pbs/test1_perl.pbs',
  'test2' => 'target/pbs/test2_perl.pbs',
};
my $actual_pbs = $test->get_pbs_files( $config, "task" );
is_deeply( $actual_pbs, $expect_pbs );

#test get_pbs_source
my $expect_pbs_source = {
  "target/pbs/test1_perl.pbs" => ["test1"],
  "target/pbs/test2_perl.pbs" => ["test2"],
};
my $actual_pbs_source = $test->get_pbs_source( $config, "task" );
is_deeply( $actual_pbs_source, $expect_pbs_source );

#test result_pbs
my $expect_result_pbs_map = {
  'test1' => 'target/pbs/test1_perl.pbs',
  'test2' => 'target/pbs/test2_perl.pbs',
};

my $result_pbs_map = $test->get_result_pbs( $config, "task" );
is_deeply( $result_pbs_map, $expect_result_pbs_map );

1;
