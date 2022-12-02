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

  initDefaultValue( $def, "MAP_EXOGENOUS", "off");
  initDefaultValue( $def, "ADAPTER_SEQ", getValue($def, "adapter", "guessKnown"));
  initDefaultValue( $def, "RANDOM_BARCODE_LENGTH", getValue($def, "fastq_remove_random", 0));
  initDefaultValue( $def, "RANDOM_BARCODE_LOCATION", "-5p -3p");
  initDefaultValue( $def, "MIN_READ_LENGTH", getValue($def, "min_read_length", 16));
  
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

  my $exceRpt_DB = getValue($def, "exceRpt_DB");
  my $host = getValue($def, "exceRpt_host");
  my $image = getValue($def, "exceRpt_image");

  my $ADAPTER_SEQ = getValue($def, "ADAPTER_SEQ");
  my $MAP_EXOGENOUS = getValue($def, "MAP_EXOGENOUS");
  my $RANDOM_BARCODE_LENGTH = getValue($def, "RANDOM_BARCODE_LENGTH");
  my $RANDOM_BARCODE_LOCATION = getValue($def, "RANDOM_BARCODE_LOCATION");
  my $MIN_READ_LENGTH = getValue($def, "MIN_READ_LENGTH");

  $config->{"exceRpt"} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/exceRpt",
    option                => "
input_dir=\$(dirname __FILE__)
input_file=\$(basename __FILE__)

echo input_dir=\$input_dir
echo input_file=\$input_file
echo sample_name=__NAME__
echo host_genome=$host

mkdir -p tmp

singularity run \\
    -c -e \\
    -B \$input_dir:/exceRptInput \\
    -B `pwd`:/exceRptOutput \\
    -B `pwd`/tmp:/tmp \\
    -B $exceRpt_DB/NCBI_taxonomy_taxdump:/exceRpt_DB/NCBI_taxonomy_taxdump \\
    -B $exceRpt_DB/$host:/exceRpt_DB/$host \\
    -B $exceRpt_DB/miRBase:/exceRpt_DB/miRBase \\
    -B $exceRpt_DB/ribosomeDatabase:/exceRpt_DB/ribosomeDatabase \\
    -B $exceRpt_DB/Genomes_BacteriaFungiMammalPlantProtistVirus:/exceRpt_DB/Genomes_BacteriaFungiMammalPlantProtistVirus \\
    $image \\
    INPUT_FILE_PATH=/exceRptInput/\$input_file \\
    MAIN_ORGANISM_GENOME_ID=$host \\
    RANDOM_BARCODE_LENGTH=$RANDOM_BARCODE_LENGTH \\
    RANDOM_BARCODE_LOCATION='$RANDOM_BARCODE_LOCATION' \\
    MIN_READ_LENGTH=$MIN_READ_LENGTH \\
    ADAPTER_SEQ=$ADAPTER_SEQ \\
    N_THREADS=8 \\
    JAVA_RAM=40G \\
    REMOVE_LARGE_INTERMEDIATE_FILES=true \\
    SAMPLE_NAME='__NAME__' \\
    MAP_EXOGENOUS=$MAP_EXOGENOUS
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
    no_output             => 1, #no output defined in command
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
