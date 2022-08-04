#!/usr/bin/perl
package CQS::ProgramIndividualWrapper;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::Task;
use File::Spec;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_piw";
  bless $self, $class;
  return $self;
}

sub getFiles {
  my ( $pfiles1, $join_arg, $first_file_only ) = @_;
  my $result = [];
  if ($join_arg) {
    push @$result, join( ',', @$pfiles1 );
  }
  elsif ($first_file_only) {
    push @$result, $pfiles1->[0];
  }
  else {
    $result = $pfiles1;
  }
  return ($result);

}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );
  my $task_suffix = $self->{_task_suffix};

  my $interpretor = get_option( $config, $section, "interpretor", "" );
  my $program = get_option( $config, $section, "program" );
  
  if (get_option($config, $section, "check_program", 1)){
    if ( !File::Spec->file_name_is_absolute($program) ) {
      $program = dirname(__FILE__) . "/$program";
    }
    if ( !( -e $program ) ) {
      die("program $program defined but not exists!");
    }
  }

  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder" );
  my $output_file_ext       = get_file_ext( $config, $section );
  my $first_file_only       = get_option( $config, $section, "first_file_only", 0 );
  my $output_arg            = get_option( $config, $section, "output_arg", "" );
  my $join_arg              = get_option( $config, $section, "join_arg", 0 );

  my ( $parameterSampleFile1, $parameterSampleFile1arg, $parameterSampleFile1JoinDelimiter ) = get_parameter_sample_files( $config, $section, "source" );
  my ( $parameterSampleFile2, $parameterSampleFile2arg, $parameterSampleFile2JoinDelimiter ) = get_parameter_sample_files( $config, $section, "parameterSampleFile2" );
  my ( $parameterSampleFile3, $parameterSampleFile3arg, $parameterSampleFile3JoinDelimiter ) = get_parameter_sample_files( $config, $section, "parameterSampleFile3" );

  $option = $option . " " . get_parameter_file_option($config, $section);

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  my $bFound2 = 0;
  my $bFound3 = 0;
  for my $sample_name ( sort keys %$parameterSampleFile1 ) {
    if (defined $parameterSampleFile2->{$sample_name}){
      $bFound2 = 1;
    }
    if (defined $parameterSampleFile3->{$sample_name}){
      $bFound3 = 1;
    }
  }
  
  my $param_option2 = "";
  if(not $bFound2){
    my $file2 = save_parameter_sample_file( $config, $section, "parameterSampleFile2", "${result_dir}/${task_name}_${task_suffix}_fileList2.list" );
    if($file2 ne ""){
      $file2 = basename($parameterSampleFile2);
      my $arg2 = get_option($config, $section, "parameterSampleFile2_arg", "");
      $param_option2 = $arg2 . " " . $file2;
    }
  }
  
  my $param_option3 = "";
  if(not $bFound3){
    my $file3 = save_parameter_sample_file( $config, $section, "parameterSampleFile3", "${result_dir}/${task_name}_${task_suffix}_fileList3.list" );
    if($file3 ne ""){
      $file3 = basename($file3);
      my $arg3 = get_option($config, $section, "parameterSampleFile3_arg", "");
      $param_option3 = $arg3 . " " . $file3;
    }
  }

  my $expect_result = $self->result($config, $section);

  for my $sample_name ( sort keys %$parameterSampleFile1 ) {
    my $pfiles1 = $parameterSampleFile1->{$sample_name};
    my $pfiles  = getFiles( $pfiles1, $join_arg, $first_file_only );
    my $idxend  = scalar(@$pfiles) - 1;

    my $cur_dir = $output_to_same_folder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $final_file = $expect_result->{$sample_name}[-1];
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    for my $i ( 0 .. $idxend ) {
      my $pfile1 = $pfiles->[$i];
      my $final_file = ( $first_file_only or ( 1 == scalar(@$pfiles) ) ) ? $sample_name . $output_file_ext : basename($pfile1) . $output_file_ext;

      my $curOption = "";
      if(not $bFound2){
        $curOption = $curOption . " " . $param_option2;
      }
      
      if(not $bFound3){
        $curOption = $curOption . " " . $param_option3;
      }
      
      if ( defined $parameterSampleFile2->{$sample_name} ) {
        $curOption = $curOption . " " . $parameterSampleFile2arg . " \"" . $parameterSampleFile2->{$sample_name}[$i] . "\"";
      }

      if ( defined $parameterSampleFile3->{$sample_name} ) {
        $curOption = $curOption . " " . $parameterSampleFile3arg . " \"" . $parameterSampleFile3->{$sample_name}[$i] . "\"";
      }

      print $pbs "
$interpretor $program $option $parameterSampleFile1arg \"$pfile1\" $curOption $output_arg $final_file

";
    }
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;
  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section, 0 );
  my $task_suffix = $self->{_task_suffix};

  my $first_file_only = get_option( $config, $section, "first_file_only", 0 );

  my $output_file  = get_option( $config, $section, "output_file", "source" );
  my ($source_files, $source_file_arg, $source_file_join_delimiter) = get_parameter_sample_files( $config, $section, $output_file );

  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder" );
  my $join_arg              = get_option( $config, $section, "join_arg", 0 );
	my $output_exts = get_output_ext_list( $config, $section );

  my $result = {};
  for my $sample_name ( sort keys %$source_files ) {
    my $cur_dir = $output_to_same_folder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );

    my $pfiles1 = $source_files->{$sample_name};
    my $pfiles  = getFiles( $pfiles1, $join_arg, $first_file_only );

    my @result_files = ();
    for my $pfile ( @pfiles ) {
      foreach my $ext (@output_exts) {
        my $other_file = $first_file_only ? $sample_name . $ext : basename($pfile) . $ext;
        push( @result_files, "${cur_dir}/$other_file" );
      }
    }

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
