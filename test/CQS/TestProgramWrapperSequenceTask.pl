#!/usr/bin/env perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use CQS::SequenceTaskSlurmSlim;
use Test::More tests => 1;

my $email = "quanhu.sheng.1\@vumc.org";

my $config = {
  general => {
    task_name => "sequencetask",
  },
  files => {
    "G1" => [ "S1_1", "S1_2" ],
  },
  "one2one1" => {
    class                 => "CQS::ProgramWrapperOneToOne",
    target_dir            => "c:/temp/one2one1",
    program               => "one2one",
    check_program         => 0,
    source_ref            => "files",
    output_arg            => "",
    output_file_prefix    => "",
    output_file_ext       => ".1.o2o1",
    output_other_ext      => ".2.o2o1",
    output_to_same_folder => 1,
    pbs                   => {
      email => $email,
    }
  },
  "one2many" => {
    class                 => "CQS::ProgramWrapperOneToMany",
    target_dir            => "c:/temp/one2many",
    source_ref            => "one2one1",
    output_arg            => "",
    output_file_prefix    => "",
    output_file_ext       => "._ITER_.1.o2m",
    output_other_ext      => "._ITER_.2.o2m",
    program               => "one2many",
    check_program         => 0,
    iteration             => 2,
    output_to_same_folder => 1,
    pbs                   => {
      email => $email,
    }
  },
  "one2one2" => {
    class                 => "CQS::ProgramWrapperOneToOne",
    target_dir            => "c:/temp/one2one2",
    program               => "one2one",
    check_program         => 0,
    source_ref            => ["one2many"],
    output_arg            => "",
    output_file_prefix    => "",
    output_file_ext       => ".1.o2o2",
    output_other_ext      => ".2.o2o2",
    output_to_same_folder => 1,
    pbs                   => {
      email => $email,
    }
  },
  "many2one" => {
    class                    => "CQS::ProgramWrapperManyToOne",
    target_dir               => "c:/temp/many2one",
    program                  => "many2one",
    check_program            => 0,
    source_ref               => ["one2one2"],
    parameterSampleFile2_ref => ["one2one1"],
    output_arg               => "",
    output_file_prefix       => ".m2o",
    output_file_ext          => ".m2o",
    output_to_same_folder    => 1,
    pbs                      => {
      email => $email,
    }
  },
  "sequencetask" => {
    class      => "CQS::SequenceTaskSlurmSlim",
    perform    => 1,
    target_dir => "c:/temp/sequencetask",
    option     => "",
    source     => {
      "step1" => [ "one2one1", "one2many", "one2one2", "many2one" ]
    },
    sh_direct => 1,
    cluster   => "slurm",
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  }
};

my $myclass = instantiate("CQS::SequenceTaskSlurmSlim");
my $deppbsmap = $myclass->get_dependent_pbs_map( $config, "sequencetask" );
is_deeply(
  $deppbsmap,
  {
    "one2one1" => { "c:/temp/one2one1/pbs/G1_o2o.pbs" => {} },
    "one2many" => {
      "c:/temp/one2many/pbs/G1_o2m.pbs" => { "c:/temp/one2one1/pbs/G1_o2o.pbs" => 1 }
    },
    "one2one2" => {
      "c:/temp/one2one2/pbs/G1_ITER_1_o2o.pbs" => { "c:/temp/one2many/pbs/G1_o2m.pbs" => 1 },
      "c:/temp/one2one2/pbs/G1_ITER_2_o2o.pbs" => { "c:/temp/one2many/pbs/G1_o2m.pbs" => 1 }
    },
    "many2one" => {
      "c:/temp/many2one/pbs/G1_m2o.pbs" => {
        "c:/temp/one2one2/pbs/G1_ITER_1_o2o.pbs" => 1,
        "c:/temp/one2one2/pbs/G1_ITER_2_o2o.pbs" => 1
      },
    },
  }
);

#performConfig($config);

1;
