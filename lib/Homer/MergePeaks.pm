#!/usr/bin/perl
package Homer::MergePeaks;

use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_mp";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  if ( defined $config->{$section}->{groups} ) {
    my %raw_files = %{ get_grouped_raw_files( $config, $section, "groups" ) };

    for my $sample_name ( sort keys %raw_files ) {
      my $pbs_file   = $self->get_pbs_filename( $pbs_dir, $sample_name );
      my $pbs_name   = basename($pbs_file);
      my $log        = $self->get_log_filename( $log_dir, $sample_name );
      my $log_desc   = $cluster->get_log_description($log);
      my $final_file = "${sample_name}.merged.peaks";

      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
      print $pbs "mergePeaks -prefix ${sample_name}.merged ";
      my @peak_files = @{ $raw_files{$sample_name} };
      for my $peak_file (@peak_files) {
        print $pbs "\\\n  \"$peak_file\" \n";
      }
      $self->close_pbs( $pbs, $pbs_file );
    }
  }
  else {
    my %raw_files = %{ get_raw_files( $config, $section ) };

    my $pbs_file   = $self->get_pbs_filename( $pbs_dir, $task_name );
    my $pbs_name   = basename($pbs_file);
    my $log        = $self->get_log_filename( $log_dir, $task_name );
    my $log_desc   = $cluster->get_log_description($log);
    my $final_file = "${task_name}.merged.peaks";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
    print $pbs "mergePeaks -prefix ${task_name}.merged";
    for my $sample_name ( sort keys %raw_files ) {
      my @peak_files = @{ $raw_files{$sample_name} };
      my $peak_file  = $peak_files[0];
      print $pbs " \"$peak_file\"";
    }

    print $pbs "\n";
    $self->close_pbs( $pbs, $pbs_file );
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  $result->{$task_name} = filter_array( ["$result_dir/${task_name}.merged.peaks"], $pattern );
  return $result;
}

1;
