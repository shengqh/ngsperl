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
  my $hasCondition = 0;
  for my $sample (sort keys %$batch_values){
    my $sample_values = $batch_values->{$sample};
    if(is_hash($sample_values)){
      for my $key (keys %$sample_values) {
        if ($key =~ "Condition") {
          $hasCondition = 1;
          last;
        }
      }
    }
  }

  if ($hasCondition){
    print $batch "Samples\tBatch\tConditions\n";
  }else{
    print $batch "Samples\tBatch\n";
  }

  for my $sample (sort keys %$batch_values){
    my $sample_values = $batch_values->{$sample};
    if(is_hash($sample_values)){
      if($hasCondition){
        print $batch "$sample\t" . $sample_values->{Batch} . "\t" . $sample_values->{Condition} . "\n";
      }else{
        print $batch "$sample\t" . $sample_values->{Batch} . "\n";
      }
    }else{
      print $batch "$sample\t" . $sample_values . "\n";
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

  my $count_task = ["mageck_count", ".count.txt"];
  if(getValue($def, "perform_batch_correction", 0)){
    defined $def->{batch} or die "define batch in configuration for batch correction.";
    my $batch_values=$def->{batch};
    my $batch_file=$target_dir . "/" . $task_name . "_batch.txt";
    writeBatchFile($batch_values, $batch_file);

    my $rm_batch_task="mageck_rm_batch";
    $config->{$rm_batch_task} = {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      docker_prefix            => "mageck_",
      rCode                    => "",
      target_dir               => "${target_dir}/$rm_batch_task",
      option                   => "",
      parameterFile1_ref => $count_task,
      parameterFile2 => $batch_file,
      rtemplate                => "../CRISPR/rm-batch.R",
      output_file              => "",
      output_file_ext          => ".corrected-count.txt",
      sh_direct                => 1,
      pbs                      => {
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    };

    push(@$summary, $rm_batch_task);
    $count_task = $rm_batch_task;
  }

  if(defined $def->{mageck_test}){
    my $cnv_option = "";
    if(getValue($def, "perform_cnv_correction", 0)){
      if (! -e $def->{cnv_norm}){
        die "define cnv_norm file for cnv correction in mageck test.";
      }
      if (! $def->{cnv_norm_cell_line}){
        die "define cnv_norm_cell_line for cnv correction in mageck test.";
      }
      $cnv_option = "--cell-line " . getValue($def, "cnv_norm_cell_line");
    }

    $config->{"mageck_test"} = {
      class                 => "CQS::ProgramWrapperOneToOne",
      perform               => 1,
      target_dir            => "$target_dir/mageck_test",
      docker_prefix         => "mageck_",
      option                => "test --pdf-report __FILE__ -n __NAME__ $cnv_option",
      interpretor           => "",
      check_program         => 0,
      program               => "mageck",
      source                => $def->{mageck_test},
      source_arg            => "",
      source_join_delimiter => "",
      output_to_same_folder => 0,
      parameterFile1_ref    => $count_task,
      parameterFile1_arg    => "-k",
      parameterFile2        => $def->{cnv_norm},
      parameterFile2_arg    => "--cnv-norm",
      output_arg            => "-n",
      output_to_folder      => 1,
      output_file_prefix    => "",
      output_file_ext       => ".count.txt",
      output_other_ext      => "",
      sh_direct             => 1,
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
