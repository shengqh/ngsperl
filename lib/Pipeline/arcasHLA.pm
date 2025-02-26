#!/usr/bin/perl
package Pipeline::arcasHLA;

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

our %EXPORT_TAGS = ( 'all' => [qw(performArcasHLA)] );

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
  my $is_paired_end = getValue($def, "is_paired_end", 1);
  my $single_end_option = $is_paired_end ? "" : "--single";
  my $output_file_ext = $is_paired_end ? ".extracted.1.fq.gz" : ".extracted.fq.gz";
  my $output_other_ext = $is_paired_end ? ".extracted.2.fq.gz" : "";
  my $min_count = getValue($def, "min_count", 75);

  $config->{"extract"} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/extract",
    init_command          => "",
    option                => "
if [[ ! -s __NAME__.bam ]]; then
  ln -s __FILE__ __NAME__.bam
fi

arcasHLA extract -t 8 -v $single_end_option --log __NAME__.log __NAME__.bam",
    docker_prefix         => "arcashla_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => "files",
    source_arg            => "",
    source_join_delimiter => "",
    output_to_same_folder => 0,
    output_arg            => "-o",
    output_to_folder      => 1,
    output_file_prefix    => "",
    output_file_ext       => $output_file_ext,
    output_other_ext      => $output_other_ext,
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };

  $config->{"genotype"} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/genotype",
    option                => "
arcasHLA genotype -t 8 -v $single_end_option --min_count $min_count --log __NAME__.log",
    docker_prefix         => "arcashla_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => "extract",
    source_arg            => "",
    source_join_delimiter => " ",
    output_to_same_folder => 1,
    output_to_folder      => 1,
    output_arg            => "-o",
    output_file_prefix    => "",
    output_file_ext       => ".genotype.json",
    output_other_ext       => ".genes.json,.alignment.p",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => "40gb"
    },
  };

  $config->{"merge"} = {
    class                 => "CQS::ProgramWrapper",
    perform               => 1,
    target_dir            => "$target_dir/merge",
    option                => "
arcasHLA merge -i $target_dir/genotype/result --run $task_name -o .

(head -n 1 __NAME__.genotypes.tsv && tail -n +2 __NAME__.genotypes.tsv | sort) > __NAME__.genotypes.sorted.tsv

mv __NAME__.genotypes.sorted.tsv __NAME__.genotypes.tsv

",
    docker_prefix         => "arcashla_",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => "genotype",
    source_arg            => "-i",
    source_join_delimiter => " ",
    output_to_same_folder => 1,
    output_to_folder      => 1,
    output_arg            => "-o",
    output_file_prefix    => "",
    output_file_ext       => ".genotypes.tsv",
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => ["extract", "genotype"],
      step2 => ["merge"],
    },
    sh_direct => 0,
    cluster   => $cluster,
    pbs       => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => $def->{sequencetask_run_time},
      "mem"       => "40gb"
    },
  };

  return ($config);
}

sub performArcasHLA {
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

sub performArcasHLATask {
  my ( $def, $task ) = @_;

  my $config = getConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
