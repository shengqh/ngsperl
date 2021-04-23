#!/usr/bin/perl
package CQS::GroupPickTask;

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
  my $groups = get_raw_files($config, $section, "groups");
  my $pick_index = get_option( $config, $section, "sample_index_in_group", 0);
  my $can_be_empty = get_option( $config, $section, "can_be_empty", 0);

  my $result = {};

  for my $group_name (keys %$groups) {
    my $samples = $groups->{$group_name};
    my $sample_name = $samples->[$pick_index];
    if($can_be_empty && (not defined $sample_name)){
      $result->{$group_name} = [];
    }else{
      my $result_files = $raw_files->{$sample_name};
      $result->{$group_name} = filter_array( $result_files, $pattern );
    }
  }

  return $result;
}


sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my $raw_files = get_raw_files($config, $section);
  my $groups = get_raw_files($config, $section, "groups");
  my $pick_index = get_option( $config, $section, "sample_index_in_group", 0);

  my $result = {};
  my $previous_section = $config->{$section}{"source_ref"};
  if (defined $config->{$previous_section}) {
    my $previous_classname = $config->{$previous_section}{class};

    if (defined $previous_classname) {
      my $myclass = instantiate($previous_classname);
      my $previous_pbs = $myclass->get_pbs_files($config, $previous_section);

      my $task_name = get_task_name( $config, $section );

      for my $group_name (keys %$groups) {
        my $samples = $groups->{$group_name};
        my $sample_name = $samples->[$pick_index];
        my $result_files = $raw_files->{$sample_name};
        if (defined $previous_pbs->{$sample_name}) {
          $result->{$group_name} = $previous_pbs->{$sample_name};
        }else{
          $result->{$group_name} = $previous_pbs->{$task_name};
        }
      }
    }
  }
  
  return $result;
}

1;
