#!/usr/bin/perl
package CQS::ProgramWrapper;

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
use File::Spec;
use Cwd;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_pw";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  my $task_suffix = get_option( $config, $section, "suffix", "" );
  $self->{_task_suffix} = $task_suffix;

  my $interpretor = get_option( $config, $section, "interpretor", "" );
  my $program = get_option( $config, $section, "program" );
  if ( !File::Spec->file_name_is_absolute($program) ) {
    $program = Cwd::realpath( dirname(__FILE__) . "/$program" );
  }
  if ( !( -e $program ) ) {
    die("program $program defined but not exists!");
  }

  my $output_ext = get_option( $config, $section, "output_ext", 0 );

  my $parameterSampleFile1 = save_parameter_sample_file( $config, $section, "parameterSampleFile1", "${result_dir}/${task_name}_${task_suffix}_fileList1.list" );
  my $parameterSampleFile2 = save_parameter_sample_file( $config, $section, "parameterSampleFile2", "${result_dir}/${task_name}_${task_suffix}_fileList2.list" );
  my $parameterSampleFile3 = save_parameter_sample_file( $config, $section, "parameterSampleFile3", "${result_dir}/${task_name}_${task_suffix}_fileList3.list" );

  my $parameterFile1 = parse_param_file( $config, $section, "parameterFile1", 0 );
  my $parameterFile2 = parse_param_file( $config, $section, "parameterFile2", 0 );
  my $parameterFile3 = parse_param_file( $config, $section, "parameterFile3", 0 );
  if ( !defined($parameterFile1) ) {
    $parameterFile1 = "";
  }
  if ( !defined($parameterFile2) ) {
    $parameterFile2 = "";
  }
  if ( !defined($parameterFile3) ) {
    $parameterFile3 = "";
  }

  my $pbs_file   = $self->get_pbs_filename( $pbs_dir, $task_name, ".pbs" );
  my $pbs_name   = basename($pbs_file);
  my $log        = $self->get_log_filename( $log_dir, $task_name, ".log" );
  my $final_file = $self->get_file( $result_dir, $task_name, $output_ext );
  my $log_desc   = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
  print $pbs "
$interpretor $program $option $parameterSampleFile1 $parameterSampleFile2 $parameterSampleFile3 $parameterFile1 $parameterFile2 $parameterFile3 $final_file
";

  $self->close_pbs( $pbs, $pbs_file );

  print "!!!pbs file $pbs_file created, you can run this pbs file to submit to cluster.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  my $task_suffix = get_option( $config, $section, "suffix", "" );
  $self->{_task_suffix} = $task_suffix;

  my $output_ext       = get_option( $config, $section, "output_ext",       "" );
  my $final_file = $self->get_file( $result_dir, $task_name, $output_ext );

  my $result = {};
  my @result_files = ();
  push( @result_files, $final_file );
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
