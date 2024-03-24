#!/usr/bin/perl
package Pipeline::CellRanger;

use strict;
use warnings;
use List::Util qw(first);
use File::Basename;
use Storable qw(dclone);
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use scRNA::Modules;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performCellRanger)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );
  initDefaultValue( $def, "perform_preprocessing",   0 );

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

  my $task_name = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  $summary = [@$individual, @$summary];

  my $target_dir = $def->{target_dir};

  $config->{count_files} = $def->{count_files};
  $config->{vdj_files} = $def->{vdj_files};

  my $count_task_name;

  if(defined $def->{multi_files} && getValue($def, "cellranger_multi_mode", 1)){
    my $multi_jobmode = getValue($def, "multi_jobmode", "local");
    $count_task_name = "multi_${multi_jobmode}";
    addCellRangerMulti($config, 
      $def, 
      $summary, 
      $target_dir, 
      $count_task_name, 
      "multi_files", 
      getValue($def, "csv_config"),
      $multi_jobmode);  
  }

  if(defined $def->{count_files}){
    my $count_jobmode = getValue($def, "count_jobmode", "local");
    $count_task_name = (defined $def->{count_chemistry}) ? "count_chemistry_$count_jobmode" : "count_$count_jobmode";
    addCellRangerCount($config, 
      $def, 
      $summary, 
      $target_dir, 
      $count_task_name, 
      getValue($def, "count_fastq_folder", ""),
      "count_files", 
      getValue($def, "count_reference"),
      $count_jobmode,
      $def->{count_chemistry});
  }

  if (defined $def->{vdj_files}) {
    my $vdj_jobmode = getValue($def, "vdj_jobmode", "local");
    my $vdj_task_name = (defined $def->{vdj_chain}) ? "vdj_chain_$vdj_jobmode" : "vdj_$vdj_jobmode";
    addCellRangerVdj($config, 
      $def, 
      $summary, 
      $target_dir, 
      $vdj_task_name, 
      getValue($def, "vdj_fastq_folder"),
      "vdj_files", 
      getValue($def, "vdj_reference"),
      $vdj_jobmode,
      $def->{vdj_chain});
  }

  my $cellranger_summary_task = "${count_task_name}_summary";
  $config->{$cellranger_summary_task} = {
    class                    => "CQS::UniqueRmd",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $cellranger_summary_task,
    report_rmd_file => "../scRNA/cellranger_summary.rmd",
    additional_rmd_files => "../CQS/reportFunctions.R",
    parameterSampleFile1_ref       => [$count_task_name, "metrics_summary.csv"],
    parameterSampleFile2_ref       => [$count_task_name, "web_summary.html"],
    output_file_ext      => ".cellranger.html",
    output_other_ext  => ".cellranger.html",
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "12",
      "mem"       => "10gb" 
    },
  };
  push( @$summary, $cellranger_summary_task );

  if(getValue($def, "perform_individual_qc", 1)){
    my $qc_pattern = "filtered_feature_bc_matrix.h5";
    #add_individual_qc($config, $def, $summary, $target_dir, $pipseeker_qc, undef, [$pipseeker, $qc_pattern], undef, undef, undef);
    my ($raw_individual_qc_task, $qc_report_task, $signacX_ref, $singleR_ref, $azimuth_ref) = add_individual_qc_tasks($config, $def, $summary, $target_dir, $task_name, "", undef, [$count_task_name, $qc_pattern], undef, undef);
  }

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => $summary,
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

sub performCellRanger {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

1;
