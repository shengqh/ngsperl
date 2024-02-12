#!/usr/bin/perl
package Pipeline::Lipidomics;

use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use List::Util qw(first);
use Storable qw(dclone);
use Hash::Merge qw( merge );

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::Preprocession;
use Pipeline::PipelineUtils;

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performLipidomics)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );

  initDefaultValue( $def, "combine_by", 'HIGHEST' ); #HIGHEST, SUM, BOTH

  initDefaultValue( $def, "max_thread",                '8' );
  initDefaultValue( $def, "sequencetask_run_time",     '12' );

  return $def;
}

sub getLipidomicsConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def->{perform_preprocessing} = 0;
  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);
  my $task_name = $def->{task_name};
  my $email = $def->{email};

  my $target_dir = $def->{target_dir};
  my $tasks = [@$individual, @$summary];

  my $preprocess_task = "preprocess";
  $config->{$preprocess_task} = {
    class => "CQS::UniqueRmd",
    target_dir => $target_dir . "/$preprocess_task",
    report_rmd_file => "../Lipidomics/Preprocess_main.Rmd",
    additional_rmd_files => "../CQS/countTableVisFunctions.R;../CQS/reportFunctions.R;../Lipidomics/lipidomics_func.R;../Lipidomics/lipidomics_preprocessing.Rmd;../Lipidomics/lipidomics_category.Rmd",
    option => "",
    parameterSampleFile1 => {
      task_name => $task_name,
      email => $email,
      affiliation => getValue($def, "affiliation", "CQS/Biostatistics, VUMC"),
      nomenclature_file => getValue($def, "nomenclature_file"),
      combine_by => getValue($def, "combine_by"),
    },
    parameterSampleFile2_ref => "files",
    suffix => ".preprocess",
    output_file_ext => ".preprocess.html,.files.txt,.pos.sample_meta.csv",
    can_result_be_empty_file => 0,
    sh_direct   => 1,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    },
  };
  push(@$tasks, $preprocess_task);

  if($def->{perform_qc}){
    my $qc_task = "qc";
    $config->{$qc_task} = {
      class => "CQS::UniqueRmd",
      target_dir => $target_dir . "/$qc_task",
      report_rmd_file => "../Lipidomics/QC_main.Rmd",
      additional_rmd_files => "../CQS/countTableVisFunctions.R;../CQS/reportFunctions.R;../Lipidomics/lipidomics_func.R",
      option => "",
      parameterSampleFile1 => {
        task_name => $task_name,
        email => $email,
        affiliation => getValue($def, "affiliation", "CQS/Biostatistics, VUMC"),
        remove_sample_pattern => $def->{remove_sample_pattern},
        remove_sample_description => $def->{remove_sample_description},
    },
      parameterSampleFile2_ref => ["preprocess", ".files.txt"],
      parameterSampleFile3_ref => ["preprocess", ".pos.sample_meta.csv"],
      suffix => ".qc",
      output_file_ext => ".qc.html,.files.txt",
      can_result_be_empty_file => 0,
      sh_direct   => 1,
      pbs => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "2",
        "mem"       => "10gb"
      },
    };
    push(@$tasks, $qc_task);

    if($def->{perform_DE_analysis}){

    }
  }

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      tasks => $tasks,
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

sub performLipidomics {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getLipidomicsConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

1;
