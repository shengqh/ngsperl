#!/usr/bin/perl
package CQS::TestProgramWrapper;

use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Hash::Merge qw( merge );
use Data::Dumper;
use Test::More tests => 2;

my $demultiplex_config = {
  general     => { task_name => "cpd" },
  demultiplex => {
    class                    => "CQS::ProgramWrapper",
    perform                  => 1,
    option                   => "demultiplex",
    target_dir               => "demultiplex",
    interpretor              => "python",
    program                  => "../CPD/analysis.py",
    parameterSampleFile1_arg => "-i",
    parameterSampleFile1     => {
      Control => ["ATCGCGAT"],
      UV      => ["GAACTGAT"]
    },
    output_arg                   => "-o",
    output_perSample_file        => "parameterSampleFile1",
    output_perSample_file_byName => 1,
    output_perSample_file_ext    => ".fastq.gz",
    output_to_same_folder        => 1,
    sh_direct                    => 1,
    pbs                          => {
      "email"     => "email",
      "emailType" => "emailType",
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  },
};

my $obj = getTaskClass( $demultiplex_config, "demultiplex" );
my $obj_result = $obj->result( $demultiplex_config, "demultiplex", "" );

is_deeply(
  $obj_result,
  {
    'Control' => [ 'demultiplex/result/Control.fastq.gz' ],
    'UV'      => [ 'demultiplex/result/UV.fastq.gz' ]
  }
);

my $pbs_result = $obj->get_pbs_files( $demultiplex_config, "demultiplex", "" );

is_deeply(
  $pbs_result,
  {
    'cpd' => 'demultiplex/pbs/cpd_pw.pbs',
  }
);

1;
