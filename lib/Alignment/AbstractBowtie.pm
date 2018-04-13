#!/usr/bin/perl
package Alignment::AbstractBowtie;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;
use Alignment::AlignmentUtils;

our @ISA = qw(CQS::Task);

sub getOutputToSameFolder {
  my ( $self, $config, $section ) = @_;
  return get_option( $config, $section, "output_to_same_folder", 1 );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $sort_by_coordinate = get_option( $config, $section, "sort_by_coordinate",    1 );
  my $mark_duplicates    = hasMarkDuplicate( $config->{$section} );
  my $outputToSameFolder = $self->getOutputToSameFolder( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my $cur_dir = $outputToSameFolder ? $result_dir : $result_dir . "/$sample_name";

    my @result_files = ();
    if ( $sort_by_coordinate && $mark_duplicates ) {
      push( @result_files, "${cur_dir}/${sample_name}.rmdup.bam" );
    }
    else {
      push( @result_files, "${cur_dir}/${sample_name}.bam" );
    }
    #push( @result_files, "${cur_dir}/${sample_name}.log" );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }

  return $result;
}

1;
