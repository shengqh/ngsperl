#!/usr/bin/perl
package Pipeline::spatialTranscriptome;

use strict;
use warnings;
use List::Util qw(first);
use File::Basename;
use Storable qw(dclone);
use File::Slurp;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use Data::Dumper;
use Hash::Merge qw( merge );
use Storable qw(dclone);
use scRNA::Modules;

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(initializeSpatialTranscriptomeDefaultOptions 
  performSpatialTranscriptome)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeSpatialTranscriptomeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );
  initDefaultValue( $def, "perform_preprocessing", 0);

  return $def;
}

sub getSpatialTranscriptome {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeSpatialTranscriptomeDefaultOptions($def);

  my $project_name = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $tasks = [@$individual, @$summary];

  my $target_dir = getValue($def, "target_dir");

  if ($def->{perform_RCTD}){
    my $rctd_task = "RCTD";

    my $binsize = "8";
    my $RCTD_thread = getValue($def, "RCTD_thread", 8);

    $config->{$rctd_task} = {
      class => "CQS::IndividualR",
      target_dir => "$target_dir/$rctd_task",
      perform => 1,
      option => "",
      rtemplate => "../scRNA/Deconvolution_functions.R,../scRNA/Deconvolution_RCTD.r",
      parameterSampleFile1_ref => "files",
      parameterSampleFile2 => {
        "bin.size" => $binsize,
        "RCTD_thread" => $RCTD_thread
      },
      parameterFile1 => getValue($def, "reference"),
      sh_direct => 1,
      no_docker => 0,
      output_ext => "_${binsize}um.rds",
      pbs => {
        "nodes"    => "1:ppn=${RCTD_thread}",
        "walltime" => "8",
        "mem"      => "40gb"
      }
    };
    push (@$tasks, $rctd_task);
  }

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => $tasks,
    },
    sh_direct => 0,
    cluster   => $cluster,
    pbs       => {
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => $def->{sequencetask_run_time},
      "mem"       => "40gb"
    },
  };

  return ($config);
}

sub performSpatialTranscriptome {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getSpatialTranscriptome($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

1;

