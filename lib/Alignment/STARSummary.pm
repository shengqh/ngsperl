#!/usr/bin/perl
package Alignment::STARSummary;

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
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_ss";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $pbs_file   = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name   = basename($pbs_file);
  my $log        = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);
  
  my $logFiles = get_raw_files( $config, $section ) ;
  
  my $mapFile = $pbs_dir . "/${task_name}_sample.list";
  writeFileList($mapFile, $logFiles, 0);
  
  my $r_script = dirname(__FILE__) . "/STARSummary.r";
  if ( !-e $r_script ) {
    die "File not found : " . $r_script;
  }
  
  my $finalFile = "${task_name}.STARSummary.tsv";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $finalFile );
  print $pbs "
R --vanilla -f $r_script --args $mapFile $finalFile
";
  $self->close_pbs( $pbs, $pbs_file );

  if ( is_linux() ) {
    chmod 0755, $pbs_file;
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $result       = {};
  my @result_files = ();
  push( @result_files, "${result_dir}/${task_name}.STARSummary.tsv" );
  push( @result_files, "${result_dir}/${task_name}.STARSummary.tsv.png" );
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
