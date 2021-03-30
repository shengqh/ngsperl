#!/usr/bin/perl
package Pipeline::CRISPRScreen;

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

our %EXPORT_TAGS = ( 'all' => [qw(
  writeBatchFile
  performCRISPRScreen
)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );
  initDefaultValue( $def, "perform_preprocessing",   1 );

  return $def;
}

sub writeBatchFile {
  my ($batch_values, $batch_file) = @_;
  open( my $batch, '>', $batch_file ) or die "Cannot create $batch_file";
  print $batch "Sample_name\tBatch\n";
  for my $batch_index (0 .. scalar(@$batch_values)){
    my $batch_samples = $batch_values->[$batch_index];
    my $idx=$batch_index+1;
    for my $sample (@$batch_samples){
      print $batch "$sample\t$idx\n";
    }
  }
  close($batch);
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

  my $task_name = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $target_dir      = $def->{target_dir};
  my $mageck_library = getValue($def, "mageck_library");


  $config->{"mageck_count"} = {
    class                 => "CQS::ProgramWrapper",
    perform               => 1,
    target_dir            => "$target_dir/mageck_count",
    docker_prefix         => "mageck_",
    #init_command          => "ln -s __FILE__ __NAME__.bam",
    option                => "count --pdf-report -n __NAME__ -l $mageck_library ",
    interpretor           => "",
    check_program         => 0,
    program               => "mageck",
    source_ref            => $source_ref,
    source_type           => "array",
    source_arg            => "--fastq",
    source_join_delimiter => " ",
    source_name_arg       => "--sample-label",
    source_name_join_delimiter => ",",
    source_name_has_comma => 1,
    output_to_same_folder => 0,
    output_arg            => "-n",
    output_to_folder      => 1,
    output_file_prefix    => "",
    output_file_ext       => ".count.txt",
    output_other_ext      => "",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };

  push(@$summary, "mageck_count");

  my $count_task = "mageck_count";
  if(defined $def->{batch}){
    my $batch_values=$def->{batch};
    my $batch_file=$target_dir . "/" . $task_name . "_batch.txt";
    writeBatchFile($batch_values, $batch_file);
  }

  # $config->{"mageck_count"} = {
  #   class                 => "CQS::ProgramWrapperOneToOne",
  #   perform               => 1,
  #   target_dir            => "$target_dir/mageck_count",
  #   docker_prefix         => "mageck_",
  #   #init_command          => "ln -s __FILE__ __NAME__.bam",
  #   option                => "count --pdf-report --sample-label \"__NAME__\" --fastq __FILE__ -n __NAME__ -l $mageck_library ",
  #   interpretor           => "",
  #   check_program         => 0,
  #   program               => "mageck",
  #   source_ref            => $source_ref,
  #   source_arg            => "",
  #   source_join_delimiter => "",
  #   output_to_same_folder => 0,
  #   output_arg            => "-n",
  #   output_to_folder      => 1,
  #   output_file_prefix    => "",
  #   output_file_ext       => ".count.txt",
  #   output_other_ext      => "",
  #   sh_direct             => 0,
  #   pbs                   => {
  #     "nodes"     => "1:ppn=8",
  #     "walltime"  => "10",
  #     "mem"       => "40gb"
  #   },
  # };

  # push(@$individual, "mageck_count");

  # $config->{"mageck_count_table"} = {
  #   class                 => "CQS::ProgramWrapper",
  #   perform               => 1,
  #   target_dir            => "$target_dir/mageck_count_table",
  #   docker_prefix         => "",
  #   option                => "",
  #   interpretor           => "python3",
  #   check_program         => 1,
  #   program               => "../CRISPR/countTable.py",
  #   source_ref            => [ "mageck_count", ".count.txt" ],
  #   source_arg            => "-i",
  #   output_arg            => "-o",
  #   output_to_folder      => 0,
  #   output_file_prefix    => ".count.txt",
  #   output_file_ext       => ".count.txt",
  #   output_other_ext      => "",
  #   sh_direct             => 0,
  #   pbs                   => {
  #     "nodes"     => "1:ppn=8",
  #     "walltime"  => "10",
  #     "mem"       => "40gb"
  #   },
  # };

  # push(@$summary, "mageck_count_table");

  if(defined $def->{mageck_test}){
    $config->{"mageck_test"} = {
      class                 => "CQS::ProgramWrapperOneToOne",
      perform               => 1,
      target_dir            => "$target_dir/mageck_test",
      docker_prefix         => "mageck_",
      option                => "test --pdf-report __FILE__ -n __NAME__ ",
      interpretor           => "",
      check_program         => 0,
      program               => "mageck",
      source                => $def->{mageck_test},
      source_arg            => "",
      source_join_delimiter => "",
      output_to_same_folder => 0,
      parameterFile1_ref    => ["mageck_count", ".count.txt"],
      parameterFile1_arg    => "-k",
      output_arg            => "-n",
      output_to_folder      => 1,
      output_file_prefix    => "",
      output_file_ext       => ".count.txt",
      output_other_ext      => "",
      sh_direct             => 0,
      pbs                   => {
        "nodes"     => "1:ppn=8",
        "walltime"  => "10",
        "mem"       => "40gb"
      },
    };

    push(@$summary, "mageck_test");
  }

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => $individual,
      step2 => $summary,
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

sub performCRISPRScreen {
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
