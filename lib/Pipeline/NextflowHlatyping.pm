#!/usr/bin/perl
package Pipeline::NextflowHlatyping;

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
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performNextflowHlatyping)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';


sub initializeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType",             "FAIL" );
  initDefaultValue( $def, "cluster",               "slurm" );
  initDefaultValue( $def, "perform_preprocessing", 0 );

  return $def;
} ## end sub initializeDefaultOptions


sub getNextflowHlatypingConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

  my $task_name = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $tasks = [ @$individual, @$summary ];

  my $target_dir = $def->{target_dir};

  my $nextflow_config  = getValue( $def, "nextflow_config" );
  my $nextflow_main_nf = getValue( $def, "nextflow_main_nf" );
  my $sh_direct        = getValue( $def, "sh_direct" );

  my $nextflow_run_mode = getValue( $def, "nextflow_run_mode", "local" );

  my $nextflow_hlatyping_task = "nextflow_hlatyping";

  if ( $nextflow_run_mode eq "slurm" ) {
    # in slurm mode, we can let nextflow to handle all samples together, the we can only need to run one main script at command line.
    $config->{$nextflow_hlatyping_task} = {
      class      => "CQS::ProgramWrapper",
      target_dir => $target_dir . "/$nextflow_hlatyping_task",
      option     => "
mkdir -p $target_dir/$nextflow_hlatyping_task/result/.mplconfig
export MPLCONFIGDIR=$target_dir/$nextflow_hlatyping_task/.mplconfig  

nextflow run $nextflow_main_nf \\
  -config $nextflow_config \\
  -profile singularity \\
  --input fileList1.list.csv \\
  --outdir . 

status=\$?
if [ \$status -ne 0 ]; then
  echo \"Error: nextflow run nf-core/latyping failed with status \$status\"
  exit \$status
else
  echo \"nextflow run nf-core/latyping completed successfully\"
  rm -rf work
fi

# for list input files: ",
      program                             => "",
      check_program                       => 0,
      parameterSampleFile1_ref            => $source_ref,
      parameterSampleFile1_header         => "sample,fastq_1,fastq_2,seq_type",
      parameterSampleFile1_join_delimiter => ",",
      parameterSampleFile1_fileFirst      => 0,
      parameterSampleFile1_col_delimiter  => ",",
      parameterSampleFile1_fileSuffix     => ".csv",
      parameterSampleFile1_suffix         => ",",
      no_prefix                           => 1,
      no_output                           => 1,
      output_ext                          => "multiqc/bismark/multiqc_report.html",
      samplename_in_result                => 0,
      output_to_same_folder               => 1,
      sh_direct                           => $sh_direct,
      no_docker                           => 1,
      pbs                                 => {
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb",
      },
    };

    push @$tasks, $nextflow_hlatyping_task;
  }

  $config->{"sequencetask"} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => { tasks => $tasks, },
    sh_direct  => 0,
    pbs        => {
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  return ($config);

} ## end sub getNextflowHlatypingConfig


sub performNextflowHlatyping {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getNextflowHlatypingConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
} ## end sub performNextflowHlatyping

1;
