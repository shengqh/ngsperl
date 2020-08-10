#!/usr/bin/perl
package CQS::ProgramWrapperOneToOne;

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
use Data::Dumper;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_o2o";
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
  my $program     = get_program( $config, $section );

  my $output_to_folder = get_option( $config, $section, "output_to_folder", 0 );
  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder", 1);
  my $output_file_prefix    = get_option( $config, $section, "output_file_prefix", (!$output_to_folder) );
  my $output_arg            = get_option( $config, $section, "output_arg" );

  my ( $parameterSampleFile1, $parameterSampleFile1arg, $parameterSampleFile1JoinDelimiter ) = get_parameter_sample_files( $config, $section, "source" );
  my ( $parameterSampleFile2, $parameterSampleFile2arg, $parameterSampleFile2JoinDelimiter ) = get_parameter_sample_files( $config, $section, "parameterSampleFile2" );
  my ( $parameterSampleFile3, $parameterSampleFile3arg, $parameterSampleFile3JoinDelimiter ) = get_parameter_sample_files( $config, $section, "parameterSampleFile3" );

  my ( $parameterFile1, $parameterFile1arg ) = get_parameter_file( $config, $section, "parameterFile1" );
  my ( $parameterFile2, $parameterFile2arg ) = get_parameter_file( $config, $section, "parameterFile2" );
  my ( $parameterFile3, $parameterFile3arg ) = get_parameter_file( $config, $section, "parameterFile3" );

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
      my $arg2 = get_option( $config, $section, "parameterSampleFile2_arg", "" );
      $param_option2 = $arg2 . " " . $file2;
    }
  }

  my $param_option3 = "";
  if ( not $bFound3 ) {
    my $file3 = save_parameter_sample_file( $config, $section, "parameterSampleFile3", "${result_dir}/${task_name}_${task_suffix}_fileList3.list" );
    if ( $file3 ne "" ) {
      $file3 = basename($file3);
      my $arg3 = get_option( $config, $section, "parameterSampleFile3_arg", "" );
      $param_option3 = $arg3 . " " . $file3;
    }
  }

  my $expect_result = $self->result( $config, $section );

  for my $sample_name ( sort keys %$parameterSampleFile1 ) {
    my $curOption = $option;

    my $cur_dir = $output_to_same_folder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $final_file = $expect_result->{$sample_name}[-1];
    my $pbs        = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    my $final_prefix = $output_to_folder ? "." : $sample_name . $output_file_prefix;

    #print Dumper($parameterSampleFile1->{$sample_name});

    if ($curOption =~ /__NAME__/){
      $curOption =~ s/__NAME__/$sample_name/g;
    }

    my $param_option1 = get_program_param( $parameterSampleFile1, $parameterSampleFile1arg, $parameterSampleFile1JoinDelimiter, $sample_name );
    if ($curOption =~ /__FILE__/){
      $curOption =~ s/__FILE__/$param_option1/g;
    }else{
      $curOption = $curOption . " " . $param_option1;
    }

    my $output_option = "$output_arg $final_prefix";
    if ($curOption =~ /__OUTPUT__/){
      $curOption =~ s/__OUTPUT__/$final_prefix/g;
      $output_option = "";
    }
    
    if ( not $bFound2 ) {
      $curOption = $curOption . " " . $param_option2;
    }
    else {
      $curOption = $curOption . " " . get_program_param( $parameterSampleFile2, $parameterSampleFile2arg, $parameterSampleFile2JoinDelimiter, $sample_name );
    }

    if ( not $bFound3 ) {
      $curOption = $curOption . " " . $param_option3;
    }
    else {
      $curOption = $curOption . " " . get_program_param( $parameterSampleFile3, $parameterSampleFile3arg, $parameterSampleFile3JoinDelimiter, $sample_name );
    }

    print $pbs "
$interpretor $program $curOption $parameterFile1arg $parameterFile1 $parameterFile2arg $parameterFile2 $parameterFile3arg $parameterFile3 $output_option

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my ( $source_files, $source_file_arg, $source_file_join_delimiter ) = get_parameter_sample_files( $config, $section, "source" );
  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder" );
  my $output_exts           = get_output_ext_list( $config, $section );

  my $result = {};
  for my $sample_name ( sort keys %$source_files ) {
    my $cur_dir = $output_to_same_folder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );

    my @result_files = ();
    for my $output_ext (@$output_exts) {
      my $result_file;
      if ( $output_ext ne "" ) {
        if ($output_ext =~ /__NAME__/) {
          $result_file = $output_ext;
          $result_file =~ s/__NAME__/$sample_name/g;
        }else{
          $result_file = $sample_name . $output_ext;
        }
        push( @result_files, "${cur_dir}/$result_file" );
      }
    }

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
