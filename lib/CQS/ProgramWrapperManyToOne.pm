#!/usr/bin/perl
package CQS::ProgramWrapperManyToOne;

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
  $self->{_suffix} = "_m2o";
  bless $self, $class;
  return $self;
}

sub my_get_joined_files {
  my ( $sample_name_map, $sample_name, $parameterSampleFile1, $join_delimiter ) = @_;
  my $pfiles                  = [];
  my $individual_sample_names = $sample_name_map->{$sample_name};
  for my $individual_sample_name (@$individual_sample_names) {
    my $p_invividual_files = $parameterSampleFile1->{$individual_sample_name};
    my $p_invividual_file  = $p_invividual_files->[0];
    push( @$pfiles, $p_invividual_file );
  }
  my $result = join( $join_delimiter, @$pfiles );
  return ($result);
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  my $task_suffix = get_option( $config, $section, "suffix", "" );
  $self->{_task_suffix} = $task_suffix;

  my $interpretor = get_option( $config, $section, "interpretor", "" );
  my $program     = get_program( $config, $section );

  my $output_to_folder = get_option( $config, $section, "output_to_folder", 0 );
  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder" );
  my $output_file_prefix    = get_option( $config, $section, "output_file_prefix" );
  my $output_arg            = get_option( $config, $section, "output_arg", "" );

  my ( $parameterSampleFile1, $parameterSampleFile1arg, $parameterSampleFile1JoinDelimiter ) = get_parameter_sample_files( $config, $section, "source" );
  my ( $parameterSampleFile2, $parameterSampleFile2arg, $parameterSampleFile2JoinDelimiter ) = get_parameter_sample_files( $config, $section, "parameterSampleFile2" );
  my ( $parameterSampleFile3, $parameterSampleFile3arg, $parameterSampleFile3JoinDelimiter ) = get_parameter_sample_files( $config, $section, "parameterSampleFile3" );

  $option = get_parameter_file_option($config, $section, $option);

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  my $bFound2 = 0;
  my $bFound3 = 0;
  for my $sample_name ( sort keys %$parameterSampleFile1 ) {
    if ( defined $parameterSampleFile2->{$sample_name} ) {
      $bFound2 = 1;
    }
    if ( defined $parameterSampleFile3->{$sample_name} ) {
      $bFound3 = 1;
    }
  }

  my $param_option2 = "";
  if ( not $bFound2 ) {
    my $file2 = save_parameter_sample_file( $config, $section, "parameterSampleFile2", "${result_dir}/${task_name}_${task_suffix}_fileList2.list" );
    if ( $file2 ne "" ) {
      $file2 = basename($parameterSampleFile2);
      $param_option2 = $parameterSampleFile2arg . " " . $file2;
    }
  }

  my $param_option3 = "";
  if ( not $bFound3 ) {
    my $file3 = save_parameter_sample_file( $config, $section, "parameterSampleFile3", "${result_dir}/${task_name}_${task_suffix}_fileList3.list" );
    if ( $file3 ne "" ) {
      $file3 = basename($file3);
      $param_option3 = $parameterSampleFile3arg . " " . $file3;
    }
  }

  my $expect_result = $self->result( $config, $section );

  my $sample_name_map = get_interation_sample_subsample_map($parameterSampleFile1);

  for my $sample_name ( sort keys %$sample_name_map ) {
    my $input_file = my_get_joined_files($sample_name_map, $sample_name, $parameterSampleFile1, $parameterSampleFile1JoinDelimiter);

    my $cur_dir = $output_to_same_folder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $final_file = $expect_result->{$sample_name}[-1];
    my $pbs        = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    my $final_prefix = $sample_name . $output_file_prefix;

    my $curOption = "";
    my $g_option = $option;

    if ( not $bFound2 ) {
      $curOption = $curOption . " " . $param_option2;
    }else{
      my $param2 = my_get_joined_files($sample_name_map, $sample_name, $parameterSampleFile2, $parameterSampleFile2JoinDelimiter);
      if ($param2 ne ""){
        if ($g_option =~ /__INPUT2__/){
          $g_option =~ s/__INPUT2__/$param2/g;
        }else{
          $curOption = $curOption . " " . $parameterSampleFile2arg . " \"" .$param2 . "\"";
        }
      }
    }

    if ( not $bFound3 ) {
      $curOption = $curOption . " " . $param_option3;
    }else{
      my $param3 = my_get_joined_files($sample_name_map, $sample_name, $parameterSampleFile3, $parameterSampleFile3JoinDelimiter);
      if ($param3 ne ""){
        if ($g_option =~ /__INPUT3__/){
          $g_option =~ s/__INPUT3__/$param3/g;
        }else{
          $curOption = $curOption . " " . $parameterSampleFile3arg . " \"" .$param3 . "\"";
        }
      }
    }

    if ($g_option =~ /__NAME__/){
      $g_option =~ s/__NAME__/$sample_name/g;
    }

    if ($g_option =~ /__INPUT__/){
      $g_option =~ s/__INPUT__/$input_file/g;
    }elsif (option_contains_arg($g_option, $parameterSampleFile1arg)) {
    }else{
      $g_option = "$g_option $parameterSampleFile1arg $input_file";
    }
    if ($g_option =~ /__OUTPUT__/){
      $g_option =~ s/__OUTPUT__/$final_prefix/g;
    }elsif (option_contains_arg($g_option, $output_arg)) {
    }else{
      $g_option = "$g_option $output_arg $final_prefix";
    }

    print $pbs "
$interpretor $program $curOption $g_option

";
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my ($source_files, $source_file_arg, $source_file_join_delimiter) = get_parameter_sample_files( $config, $section, "source" );
  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder" );
  my $output_exts = get_output_ext_list( $config, $section );

  my $result = {};
  
  my $sample_name_map = get_interation_sample_subsample_map($source_files);

  for my $sample_name ( sort keys %$sample_name_map ) {
    my $cur_dir = $output_to_same_folder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );
    my @result_files = ();
    for my $output_ext (@$output_exts) {
      if ( $output_ext ne "" ) {
        my $result_file = $sample_name . $output_ext;
        push( @result_files, "${cur_dir}/$result_file" );
      }
    }

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  #print  "task_name = " . $task_name . "\n";

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $result = {};

  my ($source_files, $source_file_arg, $source_file_join_delimiter) = get_parameter_sample_files( $config, $section, "source" );
  my $sample_name_map = get_interation_sample_subsample_map($source_files);

  for my $sample_name ( sort keys %$sample_name_map ) {
    if ( $self->acceptSample( $config, $section, $sample_name ) ) {
      $result->{$sample_name} = $self->get_pbs_filename( $pbs_dir, $sample_name );
    }
  }

  return $result;
}

sub get_pbs_source {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my ( $source_files, $source_file_arg, $source_file_join_delimiter ) = get_parameter_sample_files( $config, $section, "source" );
  my $sample_subsample_map = get_interation_sample_subsample_map($source_files);

  my $result = {};

  for my $sample_name (keys %$sample_subsample_map) {
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    $result->{$pbs_file} = $sample_subsample_map->{$sample_name};
  }

  return $result;
}

1;
