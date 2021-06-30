#!/usr/bin/perl
package CQS::SamplePickTask;

use strict;
use warnings;
use CQS::Task;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::ClassFactory;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  bless $self, $class;
  return $self;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my $raw_files = get_raw_files($config, $section);
  my $sample_name_hash = {};
  if(defined $config->{$section}{sample_names}){
    my $sample_names = get_option($config, $section, "sample_names");
    for my $sample_name (@$sample_names) {
      $sample_name_hash->{$sample_name} = 1;
    }
  }elsif(defined $config->{$section}{not_sample_names}){
    my $not_sample_names = get_option($config, $section, "not_sample_names");
    my $not_sample_name_hash = {};
    for my $sample_name (@$not_sample_names){
      $not_sample_name_hash->{$sample_name} = 1;
    }
    for my $sample_name (keys %$raw_files) {
      if (not $not_sample_name_hash->{$sample_name}){
        $sample_name_hash->{$sample_name} = 1;
      }
    }
  }else{
    die "Define sample_names or not_sample_names in section $section";
  }

  my $result = {};

  for my $sample_name (keys %$raw_files) {
    if ( $sample_name_hash->{$sample_name} ){
      my $result_files = $raw_files->{$sample_name};
      $result->{$sample_name} = filter_array( $result_files, $pattern );
    }
  }

  return $result;
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my $myresult = $self->result($config, $section);

  my $result = {};
  my $previous_section = $config->{$section}{"source_ref"};
  if (defined $config->{$previous_section}) {
    my $previous_classname = $config->{$previous_section}{class};

    if (defined $previous_classname) {
      my $myclass = instantiate($previous_classname);
      my $previous_pbs = $myclass->get_pbs_files($config, $previous_section);

      my $task_name = get_task_name( $config, $section );

      for my $sample_name (keys %$myresult) {
        $result->{$sample_name} = $previous_pbs->{$sample_name};
      }
    }
  }
  
  return $result;
}

1;
