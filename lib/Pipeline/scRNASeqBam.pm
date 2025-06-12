#!/usr/bin/perl
package Pipeline::scRNASeqBam;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performScRNASeqBam)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub getScRNASeqBamConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $task = $def->{task_name};
  my $email = $def->{email};
  my $target_dir = $def->{target_dir};
  my $cluster = getValue($def, "cluster", "slurm");

  my $config = {
    general => {
      task_name => $task,
      email => $email,
      cluster => $cluster,
    },
    files => $def->{files},
    "T01_sample_celltype_bcd" => {
      class =>"CQS::IndividualR",
      perform => 1,
      target_dir => $target_dir . "/T01_sample_celltype_bcd",
      rtemplate => "../scRNA/scRNA_func.r,../CQS/reportFunctions.R,../scRNA/sample_celltype_bcd.r",
      parameterFile1 => getValue($def, "seurat_meta_rds"),
      parameterSampleFile1_ref => "files",
      parameterSampleFile2 => {
        group_by_column => getValue($def, "group_by_column"),
        sample_column => getValue($def, "sample_column"),
        original_cell_column => getValue($def, "original_cell_column"),
      },
      output_file_ext => ".bcd.txt",
      sh_direct => 1,
      pbs => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "4",
        "mem"       => "10",
      },
    },
    "T02_sinto_filterbarcodes" => {
      class =>"CQS::ProgramWrapperOneToOne",
      perform => 1,
      target_dir => $target_dir . "/T02_sinto_filterbarcodes",
      program => "",
      check_program => 0,
      option => "
  echo 'Running Sinto filterbarcodes for __OUTPUT__ ...';
  sinto filterbarcodes --nproc 8 -b __FILE__ -c __FILE2__ --outdir .

  ",
      parameterSampleFile1_ref => "files",
      parameterSampleFile2_ref => "T01_sample_celltype_bcd",
      output_arg => "",
      output_file => "parameterSampleFile1",
      output_file_prefix => "",
      output_file_ext => ".bam",
      output_to_same_folder => 0,
      sh_direct => 1,
      pbs => {
        "nodes"     => "1:ppn=8",
        "walltime"  => "4",
        "mem"       => "10",
      },
    },
    sequencetask => {
      class      => getSequenceTaskClassname($cluster),
      perform    => 1,
      target_dir => "${target_dir}/sequencetask",
      option     => "",
      source     => {
        step1 => [ "T01_sample_celltype_bcd", "T02_sinto_filterbarcodes" ],
      },
      sh_direct => 0,
      pbs       => {
        "nodes" => "1:ppn=" . getValue($def, "max_thread", 8),
        "walltime" => getValue($def, "sequencetask_run_time", 48),
        "mem" => "40gb",
      },
    },
  };

  return( $config );
}

sub performScRNASeqBam {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getScRNASeqBamConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

1;
