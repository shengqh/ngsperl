#!/usr/bin/perl
package SmallRNA::miRge;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_miRge";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $mirge = get_param_file( $config->{$section}{"miRge_script"}, "miRge_script", 1, not $self->using_docker() );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  
  my @raw_file_list;
  for my $value (values %raw_files){
    push(@raw_file_list, $value->[0])
  }
  my $fastqFiles = join(',', @raw_file_list);

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $outputfile = $result_dir . "/miR.Counts.csv";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $outputfile );
  print $pbs "perl $mirge $option --CPU $thread --output . --SampleFiles $fastqFiles";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $result = {};

  my @result_files = ();
  my $outputfile = $result_dir . "/miR.Counts.csv";
  push( @result_files, $outputfile );
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
