#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use Test::More tests => 4;

my $rootDir = "C:/temp";
my $config  = {
  general => {
    task_name => "test",
    cluster   => "slurm",
  },
  "perl1" => {
    class      => "CQS::Perl",
    perform    => 1,
    target_dir => "$rootDir/perl1",
    option     => "",
    perlFile   => "../SmallRNA/filterTrnaXml.py",
    source     => {
      "test1" => [ "$rootDir/test1_1.zip", "$rootDir/test1_2.zip" ],
      "test2" => [ "$rootDir/test2_1.zip", "$rootDir/test2_2.zip" ],
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
    target_dir            => "$rootDir/perl2",
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
    class                 => "CQS::UniqueR",
    perform               => 1,
    target_dir            => "$rootDir/r1",
    option                => "",
    rtemplate             => "../SmallRNA/filterTrnaXml.py",
    source_ref            => ["perl2"],
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
  sequencetask => {
    class      => "CQS::SequenceTaskSlurmSlim",
    perform    => 1,
    target_dir => "$rootDir/sequencetask",
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

my $myclass    = instantiate( $config->{perl1}{class} );
my $expect_pbs = {
  "C:/temp/perl1/pbs/test1_perl.pbs" => ["test1"],
  "C:/temp/perl1/pbs/test2_perl.pbs" => ["test2"]
};
my $actual_pbs = $myclass->get_pbs_source( $config, "perl1" );
is_deeply( $actual_pbs, $expect_pbs );

$myclass    = instantiate( $config->{r1}{class} );
$expect_pbs = { "C:/temp/r1/pbs/test_uniqueR.pbs" => [ "test", "test1", "test2" ], };
$actual_pbs = $myclass->get_pbs_source( $config, "r1" );
is_deeply( $actual_pbs, $expect_pbs );

1;
