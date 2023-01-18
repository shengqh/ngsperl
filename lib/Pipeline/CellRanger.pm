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

  if(defined $def->{count_files}){
    if(getValue($def, "cellranger_multi_mode", 0)){
      my $multi_jobmode = getValue($def, "multi_jobmode", "local");
      my $multi_task_name = "multi_${multi_jobmode}";
      addCellRangerMulti($config, 
        $def, 
        $summary, 
        $target_dir, 
        $multi_task_name, 
        getValue($def, "count_fastq_folder"),
        "count_files", 
        getValue($def, "csv_config"),
        $multi_jobmode);
    }else{
      my $count_jobmode = getValue($def, "count_jobmode", "local");
      my $count_task_name = (defined $def->{count_chemistry}) ? "count_chemistry_$count_jobmode" : "count_$count_jobmode";
      addCellRangerCount($config, 
        $def, 
        $summary, 
        $target_dir, 
        $count_task_name, 
        getValue($def, "count_fastq_folder"),
        "count_files", 
        getValue($def, "count_reference"),
        $count_jobmode,
        $def->{count_chemistry});
    }
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
