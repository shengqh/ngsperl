#!/usr/bin/perl
package GATK::GATK3Utils;

use strict;
use warnings;
use File::Basename;
use List::Util qw[min];
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw( addCombineVariants )] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

sub addCombineVariants {
  my ($config, $def, $tasks, $target_dir, $task_name, $source_ref) = @_;

  $config->{ $task_name } = {
    class       => "CQS::ProgramWrapperOneToMany",
    perform     => 1,
    target_dir  => "$target_dir/$task_name",
    option      => "",
    interpretor => "python3",
    program     => "../Format/splitFastq.py",
    source_ref      => $source_ref,
    source_arg            => "-i",
    source_join_delimiter => ",",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_prefix    => "",
    output_file_ext       => "._ITER_.1.fastq.gz",
    output_other_ext      => "._ITER_.2.fastq.gz",
    iteration_arg         => "--trunk",
    iteration             => getValue($def, "aligner_scatter_count"),
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };

  push @$tasks, (  $split_fastq );
}
