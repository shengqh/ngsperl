#!/usr/bin/perl
package SmallRNA::SmallRNAEndPosition;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::UniqueTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_sp";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $raw_files = get_raw_files( $config, $section );

  my $file_list_file = $self->get_file( $result_dir, $task_name, ".filelist" );
  writeFileList( $file_list_file, $raw_files, 0 );

  my $py_script = dirname(__FILE__) . "/smallRNAEndPosition.py";
  if ( !-e $py_script ) {
    die "File not found : " . $py_script;
  }

  my $final_file = $task_name . ".endpoint.txt";

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
  print $pbs "python3 $py_script -i $file_list_file -o $final_file \n";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $result       = {};
  my $final_file   = $task_name . ".endpoint.txt";
  my @result_files = ($final_file);
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
