#!/usr/bin/perl
package Encode::FindPeakTask;

use strict;
use warnings;
use CQS::Task;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::ClassFactory;
use Data::Dumper;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  bless $self, $class;
  return $self;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my $raw_folders = get_raw_files($config, $section);
  my $replicates = get_raw_files($config, $section, "replicates");
  #print(Dumper($replicates) . "\n");

  my $result = {};

  for my $rep_name (keys %$replicates) {
    my $samples = $replicates->{$rep_name};
    my $rep_folder = $raw_folders->{$rep_name}[0];

    my $idx = 0;
    for my $sample_name (@$samples){
      $idx = $idx + 1;
      my $bam_files = [$rep_folder . "peak/rep" . $idx . "/" . $sample_name . "_clipped.1.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz"];
      $result->{$sample_name} = filter_array( $bam_files, $pattern );
    }
  }

  return $result;
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my $raw_folders = get_raw_files($config, $section);
  my $replicates = get_raw_files($config, $section, "replicates");

  my $result = {};
  my $previous_section = $config->{$section}{"source_ref"};
  if (defined $config->{$previous_section}) {
    my $previous_classname = $config->{$previous_section}{class};

    if (defined $previous_classname) {
      my $myclass = instantiate($previous_classname);
      my $previous_pbs = $myclass->get_pbs_files($config, $previous_section);
      #print(Dumper($previous_pbs));

      my $task_name = get_task_name( $config, $section );

      for my $rep_name (keys %$replicates) {
        my $samples = $replicates->{$rep_name};
        for my $sample_name (@$samples){
          defined $previous_pbs->{$rep_name} or die "Cannot find PBS for " . $rep_name;
          $result->{$sample_name} = $previous_pbs->{$rep_name};
        }
      }
    }
  }
  
  return $result;
}

1;
