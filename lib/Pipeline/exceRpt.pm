#!/usr/bin/perl
package Pipeline::exceRpt;

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

our %EXPORT_TAGS = ( 'all' => [qw(performExceRpt)] );

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

  my $target_dir      = $def->{target_dir};

  my $host = getValue($def, "exceRpt_host");
  my $image = getValue($def, "exceRpt_image");
  my $fastq_remove_random = getValue($def, "fastq_remove_random");
  my $min_read_length = getValue($def, "min_read_length", 16);
  my $MAP_EXOGENOUS =  getValue($def, "MAP_EXOGENOUS", "miRNA");
  
  $config->{"exceRpt"} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/exceRpt",
    option                => "
input_dir=\$(dirname __INPUT__)
input_name=\$(basename __INPUT__)

singularity run \\
    -B \$input_dir:/exceRptInput \
    -B `pwd`:/exceRptOutput \\
    -B /scratch/h_vangard_1/references/exceRpt/$host:/exceRpt_DB/$host \\
    $image \\
    INPUT_FILE_PATH=/exceRptInput/\$input_name \\
    MAIN_ORGANISM_GENOME_ID=$host \\
    RANDOM_BARCODE_LENGTH=$fastq_remove_random \\
    MIN_READ_LENGTH=$min_read_length \\
    N_THREADS=8 \\
    SAMPLE_NAME=__NAME__ \\
    MAP_EXOGENOUS=miRNA 
",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => "files",
    source_arg            => "",
    source_join_delimiter => "",
    output_to_same_folder => 0,
    output_arg            => "",
    output_file_prefix    => "",
    output_file_ext       => ".extracted.1.fq.gz",
    output_other_ext      => ".extracted.2.fq.gz",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => "40gb"
    },
  };

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      tasks => ["exceRpt"],
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

sub performExceRpt {
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

sub performExceRptTask {
  my ( $def, $task ) = @_;

  my $config = getConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
