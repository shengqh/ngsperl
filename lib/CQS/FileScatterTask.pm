#!/usr/bin/perl
package CQS::FileScatterTask;

use strict;
use warnings;
use CQS::Task;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  bless $self, $class;
  return $self;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my $iteration_zerobased = get_option( $config, $section, "iteration_zerobased", 0 );
  my $max_length = get_option( $config, $section, "max_length", 2 );
  my $step = get_option( $config, $section, "step", 1 );

  my $raw_files = get_raw_files( $config, $section, "source" );

  my $iter_start = $iteration_zerobased ? 0: 1;
  my $result = {};
  for my $sample_name ( sort keys %$raw_files ) {
    my $files = $raw_files->{$sample_name};
    my $index = 0;
    my $iter = $iter_start;
    while($index < scalar(@$files)){
      my $key = $sample_name . "_ITER_" . left_pad($iter, $max_length);
      my @result_files = ();
      my $s = 0;
      while ($s < $step){
        my $cur_index = $index + $s;
        my $result_file = $files->[$cur_index];
        push (@result_files, $result_file);
        $s = $s + 1;
      }
      $index = $index + $step;
      $result->{$key} = filter_array( \@result_files, $pattern );
      $iter = $iter + 1;
    }
  }
  return $result;
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my $raw_files = get_raw_files($config, $section);

  my $result = {};

  #should trace the previous pbs, let's hold it now.
  
  return $result;
}

1;
