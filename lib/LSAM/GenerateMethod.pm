#!/usr/bin/perl
package LSAM::GenerateMethod;

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

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_gm";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $templates = get_raw_files( $config, $section );
  my $originalName = get_option( $config, $section, "original_name" );
  my $targetNames = get_option( $config, $section, "target_names" );

  my $pbs_file = $self->get_file( $pbs_dir, $task_name, ".bat" );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $py_script = dirname(__FILE__) . "/generateMethod.py";
  if ( !-e $py_script ) {
    die "File not found : " . $py_script;
  }

  my $pbs = $self->open_pbs( $pbs_file, "", "", $path_file, $result_dir );

  for my $name (keys %$templates ) {
    my $template = $templates->{$name}[0];
    for my $targetName (@$targetNames){
      my $finalFile = $name . "_" . $targetName . ".txt";
      
      print $pbs "
python3 $py_script -i $template -o $finalFile --originalName $originalName --targetName $targetName 
";  
    }
  }
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section, 0 );

  my $templates = get_raw_files( $config, $section );
  my $targetNames = get_option( $config, $section, "target_names" );

  $result_dir =~ s/\//\\/g;
  my $result = {};
  for my $name (keys %$templates ) {
    my $template = $templates->{$name}[0];
    for my $targetName (@$targetNames){
      my $finalFile = $name . "_" . $targetName . ".txt";
      my @result_files = ();
      push( @result_files, $result_dir . "\\" . $finalFile );
      $result->{$name . "_" . $targetName} = filter_array( \@result_files, $pattern );
    }
  }

  return $result;
}
1;
