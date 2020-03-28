#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::UniqueTask;
use Test::More tests => 5;

my $config  = {
  general => {
    task_name => "test",
    cluster   => "slurm",
  },
  "perl1" => {
    class      => "CQS::Perl",
    perform    => 1,
    target_dir => "perl1",
    option     => "",
    perlFile   => "../SmallRNA/filterTrnaXml.py",
    source     => {
      "test1" => [ "test1_1.zip", "test1_2.zip" ],
      "test2" => [ "test2_1.zip", "test2_2.zip" ],
    },
    source_arg            => "-i",
    source_join_delimiter => ",",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_prefix    => "",
    output_file_ext       => "._ITER_.1.fastq.gz",
    output_other_ext      => "._ITER_.2.fastq.gz",
    iteration             => 3,
    sh_direct             => 1,
    pbs                   => {
      "email"     => "quanhu.sheng.1\@vumc.org",
      "emailType" => "FAIL",
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  },
  "perl2" => {
    class                 => "CQS::Perl",
    perform               => 1,
    target_dir            => "perl2",
    option                => "",
    perlFile              => "../SmallRNA/filterTrnaXml.py",
    source_ref            => ["perl1"],
    source_arg            => "-i",
    source_join_delimiter => ",",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_prefix    => "",
    output_file_ext       => "._ITER_.1.fastq.gz",
    output_other_ext      => "._ITER_.2.fastq.gz",
    iteration             => 3,
    sh_direct             => 1,
    pbs                   => {
      "email"     => "quanhu.sheng.1\@vumc.org",
      "emailType" => "FAIL",
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  },
  "r1" => {
    class                    => "CQS::UniqueTask",
    perform                  => 1,
    target_dir               => "r1",
    option                   => "",
    rtemplate                => "../SmallRNA/filterTrnaXml.py",
    source_ref               => ["perl2"],
    parameterSampleFile2_ref => ["perl1"],
    output_file_ext          => ".png;.csv",
    sh_direct                => 1,
    pbs                      => {
      "email"     => "quanhu.sheng.1\@vumc.org",
      "emailType" => "FAIL",
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTaskSlurmSlim",
    perform    => 1,
    target_dir => "sequencetask",
    option     => "",
    source     => {

      #step1 => ["perl1", "perl2", "r1"],
      step1 => [ "perl1", "r1" ],
    },
    sh_direct => 0,
    cluster   => "slurm",
    pbs       => {
      "email"     => "quanhu.sheng.1\@vumc.org",
      "emailType" => "FAIL",
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  }
};

my $test = CQS::UniqueTask->new();

#test result
eval {
  $test->result( $config, "r1" ); 
};
ok($@, 'exception has been thrown');
like($@, qr/Override get_result_files/);

#test pbs
my $expect_pbs = {
  'test' => 'r1/pbs/test.pbs',
};
my $actual_pbs = $test->get_pbs_files( $config, "r1" );
is_deeply( $actual_pbs, $expect_pbs );

#test get_pbs_source
my $expect_pbs_source = { "r1/pbs/test.pbs" => [ "test", "test1", "test2" ], };
my $actual_pbs_source = $test->get_pbs_source( $config, "r1" );
is_deeply( $actual_pbs_source, $expect_pbs_source );

#test result_pbs
my $expect_result_pbs_map = {
  'test' => 'r1/pbs/test.pbs',
};

my $result_pbs_map = $test->get_result_pbs($config, "r1");
is_deeply( $result_pbs_map, $expect_result_pbs_map );

1;
