#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $config = {
  general => {
    task_name => "TestOneToMany",
  },
  "OneToMany" => {
    class       => "CQS::ProgramWrapperOneToMany",
    perform     => 1,
    target_dir  => "C:/temp/1_OneToMany",
    option      => "",
    interpretor => "python",
    program     => "../SmallRNA/filterTrnaXml.py",
    source      => {
      "test1" => [ "C:/temp/test1_1.zip", "C:/temp/test1_2.zip" ],
      "test2" => [ "C:/temp/test2_1.zip", "C:/temp/test2_2.zip" ],
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
  "OneToOne" => {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "C:/temp/2_OneToOne",
    option                => "",
    interpretor           => "",
    check_program         => 0,
    program               => "bwa",
    source_ref            => "OneToMany",
    source_arg            => "-i",
    source_join_delimiter => ",",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_prefix    => ".bam",
    output_file_ext       => ".bam",
    #is_source_one2many    => 1,
    sh_direct             => 1,
    pbs                   => {
      "email"     => "quanhu.sheng.1\@vumc.org",
      "emailType" => "FAIL",
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  },
  "ManyToOne" => {
    class                 => "CQS::ProgramWrapperManyToOne",
    perform               => 1,
    target_dir            => "C:/temp/3_ManyToOne",
    option                => "",
    interpretor           => "python",
    program               => "../SmallRNA/filterTrnaXml.py",
    source_ref            => "OneToOne",
    source_arg            => "-i",
    source_join_delimiter => ",",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_prefix    => ".bam",
    output_file_ext       => ".bam",
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
    target_dir => "C:/temp/sequencetask",
    option     => "",
    source     => {
      step1 => ["OneToMany"],
      step2 => ["OneToOne"],
      step3 => ["ManyToOne"],
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

performTask($config, "sequencetask");

1;
