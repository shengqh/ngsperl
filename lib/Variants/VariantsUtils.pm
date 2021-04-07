#!/usr/bin/perl
package Variants::VariantsUtils;

use strict;
use warnings;
use File::Basename;
use List::Util qw[min];
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;
use Pipeline::PipelineUtils;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw( 
  add_combine_mutect
  add_post_mutect )] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

sub add_combine_mutect {
  my ($config, $def, $tasks, $target_dir, $task_name, $source_ref) = @_;

  $config->{$task_name} = {
    class                 => "CQS::ProgramWrapper",
    perform               => 1,
    target_dir            => "${target_dir}/${task_name}",
    option                => "",
    interpretor           => "python3",
    program               => "../GATK/mergeMutect.py",
    check_program         => 1,
    parameterSampleFile1_arg    => "-i",
    parameterSampleFile1_ref    => $source_ref,
    parameterSampleFile1_fileonly  => 0,
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_ext       => "_pass.combined.vcf",
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };

  push @$tasks, (  $task_name );
}

sub add_post_mutect {
  my ($config, $def, $target_dir, $tasks, $mutect_prefix, $mutect_index_dic, $mutect_index_key, $mutect_ref ) = @_;

  my $combineVariantsName = $mutect_prefix . getNextIndex($mutect_index_dic, $mutect_index_key) . "_merge";
  add_combine_mutect($config, $def, $tasks, $target_dir, $combineVariantsName, $mutect_ref);

  my $filterVariantsName = $mutect_prefix . getNextIndex($mutect_index_dic, $mutect_index_key) . "_filterDepth";
  $config->{$filterVariantsName} = {
    class                 => "CQS::ProgramWrapper",
    perform               => 1,
    target_dir            => "${target_dir}/${filterVariantsName}",
    option                => "",
    interpretor           => "python3",
    program               => "../GATK/filterMutect.py",
    check_program         => 1,
    parameterFile1_arg    => "-i",
    parameterFile1_ref    => [ $combineVariantsName ],
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_ext       => ".filtered.vcf",
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };
  push @$tasks, $filterVariantsName;

  if ( $def->{perform_annovar} ) {
    my $annovar_name = addAnnovar( $config, $def, $tasks, $target_dir, $filterVariantsName, ".vcf\$", $mutect_prefix, $mutect_index_dic, $mutect_index_key );
    if ( $def->{annovar_param} =~ /exac/ || $def->{annovar_param} =~ /1000g/ || $def->{annovar_param} =~ /gnomad/ ) {
      my $annovar_filter_name = addAnnovarFilter( $config, $def, $tasks, $target_dir, $annovar_name, $mutect_prefix, $mutect_index_dic, $mutect_index_key);

      if ( defined $def->{annotation_genes} ) {
        addAnnovarFilterGeneannotation( $config, $def, $tasks, $target_dir, $annovar_filter_name );
      }

      my ($annovarMaf,$annovarMafReport)=addAnnovarMafReport($config, $def, $tasks, $target_dir, $annovar_filter_name, $mutect_prefix, $mutect_index_dic, $mutect_index_key);
      return($annovarMaf,$annovarMafReport)
    }
  }
}