#!/usr/bin/perl
package CQS::SequenceTaskNoDependence;

use strict;
use warnings;
use File::Basename;
use File::Copy;
use Data::Dumper;
use CQS::Task;
use CQS::PBS;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_st";
  $self->{_forbid_tmp_folder} = 1;
  bless $self, $class;
  return $self;
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my $result = [];

  my %step_map = %{ get_raw_files( $config, $section ) };

  for my $step_name ( sort keys %step_map ) {
    my @tasks = @{ $step_map{$step_name} };

    my $samples = {};
    my $taskpbs = {};
    for my $task_section (@tasks) {
      my $classname = $config->{$task_section}{class};
      if ( !defined $classname ) {
        die "$task_section is not a valid task section.";
      }
      my $myclass = instantiate($classname);
      my $pbs_file_map = $myclass->get_pbs_files( $config, $task_section );
      for my $sample ( sort keys %{$pbs_file_map} ) {
        my $sample_files = $pbs_file_map->{$sample};
        push(@$result, $sample_files);
      }
    }
  }
  return $result;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section );

  my $cluster = get_cluster( $config, $section );

  my $pbs_files = $self->get_pbs_files($config, $section);

  my $final_name     = $task_name . "_pipeline";
  my $final_pbs      = $self->get_pbs_filename( $pbs_dir, $final_name );
  my $final_log      = $self->get_log_filename( $log_dir, $final_name );
  my $final_log_desp = $cluster->get_log_description($final_log);

  my $final = $self->open_pbs( $final_pbs, $pbs_desc, $final_log_desp, $path_file, $pbs_dir );
  for my $pbs_file (@$pbs_files){
    $final->print("sh " . $pbs_file . "\n");
  }
  $self->close_pbs( $final, $final_pbs );
}

1;
