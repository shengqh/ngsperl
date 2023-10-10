#!/usr/bin/perl
package Pipeline::EncodeATACseq;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use Pipeline::WdlPipeline;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performEncodeATACseq performEncodeATACseqTask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  initDefaultValue( $def, "cluster",               "slurm" );
  initDefaultValue( $def, "max_thread",            8 );
  initDefaultValue( $def, "sequencetask_run_time", 12 );

  #Tasks
  initDefaultValue( $def, "sra_to_fastq", 0 );

  initDefaultValue( $def, "fastq_remove_N", 0 );

  initDefaultValue( $def, "perform_cutadapt",    0 );
  if ( getValue( $def, "perform_cutadapt" ) ) {
    initDefaultValue( $def, "adapter",         "CTGTCTCTTATACACATCT" );
    initDefaultValue( $def, "min_read_length", 36 );
    initDefaultValue( $def, "cutadapt_option", "-q 30" );
  }

  if ( !defined $def->{"treatments"} ) {
    my $files = getValue( $def, "files" );
    my $groups = {};
    for my $sample_name ( sort keys %$files ) {
      $groups->{$sample_name} = [$sample_name];
    }
    $def->{"treatments"} = $groups;
  }

  initDefaultValue( $def, "perform_croo_qc", 0 );

  initDefaultValue( $def, "perform_report", 0 );
  initDefaultValue( $def, "encode_option", "" );
  #$def->{encode_option} = "-b slurm"; #or empty to use local

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  $def->{replicates} = $config->{groups};

  $def->{adapter} = "";
  
  addEncodeATACseq($config, $def, $individual, $target_dir, $source_ref, "encode_atacseq");

  addSequenceTask($config, $def, $individual, $target_dir, $summary);

  return($config);
}

sub performEncodeATACseq {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getConfig($def);
  #print(Dumper($def));

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }
  return $config;
}

sub performEncodeATACseqTask {
  my ( $def, $task ) = @_;
  my $config = getConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
