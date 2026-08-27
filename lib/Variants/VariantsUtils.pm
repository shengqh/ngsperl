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

our %EXPORT_TAGS = (
  'all' => [
    qw(
        add_combine_mutect
        add_post_mutect )
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );


sub add_combine_mutect {
  my ( $config, $def, $tasks, $target_dir, $task_name, $source_ref ) = @_;

  $config->{$task_name} = {
    class                         => "CQS::ProgramWrapper",
    perform                       => 1,
    target_dir                    => "${target_dir}/${task_name}",
    option                        => "",
    interpretor                   => "python3",
    program                       => "../GATK/mergeMutect.py",
    check_program                 => 1,
    parameterSampleFile1_arg      => "-i",
    parameterSampleFile1_ref      => $source_ref,
    parameterSampleFile1_fileonly => 0,
    output_to_same_folder         => 1,
    output_arg                    => "-o",
    output_file_ext               => "_pass.combined.vcf",
    sh_direct                     => 1,
    pbs                           => {
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  };

  push @$tasks, ($task_name);
} ## end sub add_combine_mutect


sub add_post_mutect {
  my ( $config, $def, $target_dir, $tasks, $mutect_prefix, $mutect_index_dic, $mutect_index_key, $mutect_ref ) = @_;

  my $filtered_task = undef;

  my $mutect2_individual_filtering = getValue( $def, "mutect2_individual_filtering", 0 );
  if ($mutect2_individual_filtering) {
    my $filter_task = $mutect_prefix . getNextIndex( $mutect_index_dic, $mutect_index_key ) . "_filterDepth";
    $config->{$filter_task} = {
      class         => "CQS::ProgramWrapperOneToOne",
      perform       => 1,
      target_dir    => "${target_dir}/${filter_task}",
      option        => "",
      interpretor   => "python",
      program       => "../GATK/filterIndividualMutect.py",
      check_program => 1,
      option        => " \\
      --filter_fisher_test True \\
      -n __NAME__ \\
      -i __FILE__ \\
      -o __NAME__.filtered.vcf 
    ",
      parameterSampleFile1_arg => "-i",
      parameterSampleFile1_ref => $mutect_ref,
      output_to_same_folder    => 1,
      output_arg               => "-o",
      output_file_prefix       => ".filtered.vcf",
      output_file_ext          => ".filtered.vcf",
      sh_direct                => 1,
      pbs                      => {
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    };
    push @$tasks, $filter_task;

    my $combine_task = $mutect_prefix . getNextIndex( $mutect_index_dic, $mutect_index_key ) . "_merge";
    add_combine_mutect( $config, $def, $tasks, $target_dir, $combine_task, $filter_task );

    $filtered_task = $combine_task;
  } ## end if ($mutect2_individual_filtering)
  else {
    my $combine_task = $mutect_prefix . getNextIndex( $mutect_index_dic, $mutect_index_key ) . "_merge";
    add_combine_mutect( $config, $def, $tasks, $target_dir, $combine_task, $mutect_ref );

    my $filter_task = $mutect_prefix . getNextIndex( $mutect_index_dic, $mutect_index_key ) . "_filterDepth";
    $config->{$filter_task} = {
      class                 => "CQS::ProgramWrapper",
      perform               => 1,
      target_dir            => "${target_dir}/${filter_task}",
      option                => "",
      interpretor           => "python3",
      program               => "../GATK/filterMutect.py",
      check_program         => 1,
      parameterFile1_arg    => "-i",
      parameterFile1_ref    => [$combine_task],
      output_to_same_folder => 1,
      output_arg            => "-o",
      output_file_ext       => ".filtered.vcf",
      sh_direct             => 1,
      pbs                   => {
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    };
    push @$tasks, $filter_task;

    $filtered_task = $filter_task;
  } ## end else [ if ($mutect2_individual_filtering)]

  if ( $def->{perform_annovar} ) {
    my $annovar_task = addAnnovar( $config, $def, $tasks, $target_dir, $filtered_task, ".vcf\$", $mutect_prefix, $mutect_index_dic, $mutect_index_key );
    if ( $def->{annovar_param} =~ /exac/ || $def->{annovar_param} =~ /1000g/ || $def->{annovar_param} =~ /gnomad/ ) {
      my $annovar_filter_task = addAnnovarFilter( $config, $def, $tasks, $target_dir, $annovar_task, $mutect_prefix, $mutect_index_dic, $mutect_index_key );

      if ( defined $def->{annotation_genes} ) {
        addAnnovarFilterGeneannotation( $config, $def, $tasks, $target_dir, $annovar_filter_task );
      }

      my ( $annovarMaf, $annovarMafReport ) = addAnnovarMafReport( $config, $def, $tasks, $target_dir, $annovar_filter_task, $mutect_prefix, $mutect_index_dic, $mutect_index_key );

      if ( $def->{perform_mutect2_qc} ) {
        my $qc_task = $mutect_prefix . getNextIndex( $mutect_index_dic, $mutect_index_key ) . "_qc";
        $config->{$qc_task} = {
          class                    => "CQS::UniqueRmd",
          target_dir               => "$target_dir/$qc_task",
          perform                  => 1,
          report_rmd_file          => "../Variants/variantcalling_qc.Rmd",
          additional_rmd_files     => "../CQS/reportFunctions.R",
          parameterSampleFile1_ref => $mutect_ref,
          parameterSampleFile2_ref => [ $filtered_task,       ".vcf\$" ],
          parameterSampleFile3_ref => [ $annovar_filter_task, ".filtered.tsv\$" ],
          output_file_ext          => ".variant_qc.html",
          output_other_ext         => ".variant_qc.html",
          pbs                      => {
            "nodes"    => "1:ppn=1",
            "walltime" => "24",
            "mem"      => "40gb"
          },
        };
        push @$tasks, $qc_task;
      } ## end if ( $def->{perform_mutect2_qc...})

      return ( $annovarMaf, $annovarMafReport );
    } ## end if ( $def->{annovar_param...})
  } ## end if ( $def->{perform_annovar...})

} ## end sub add_post_mutect
