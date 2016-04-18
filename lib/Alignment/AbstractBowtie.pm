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

our @ISA = qw(CQS::Task);

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $samformat = get_option( $config, $section, "samformat", 1 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my $cur_dir = $result_dir . "/$sample_name";

    my $final_file;
    if ($samformat) {
      $final_file = $sample_name . ".bam";
    }
    else {
      $final_file = $sample_name . ".out";
    }
    my @result_files = ();
    push( @result_files, $cur_dir . "/" . $final_file );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }

  return $result;
}

1;
