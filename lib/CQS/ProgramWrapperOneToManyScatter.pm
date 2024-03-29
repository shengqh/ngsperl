#!/usr/bin/perl
package CQS::ProgramWrapperOneToManyScatter;

use strict;
use warnings;
use File::Basename;
use Data::Dumper;
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
  $self->{_suffix} = "_o2m";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  my $task_suffix = get_option( $config, $section, "suffix", "" );
  $self->{_task_suffix} = $task_suffix;

  my $interpretor = get_option( $config, $section, "interpretor", "" );
  my $program     = get_program( $config, $section );

  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder" );
  my $output_file_prefix    = get_option( $config, $section, "output_file_prefix" );
  my $output_arg      = get_option( $config, $section, "output_arg" );

  my ( $parameterSampleFile1, $parameterSampleFile1arg, $parameterSampleFile1JoinDelimiter ) = get_parameter_sample_files( $config, $section, "source" );
  my ( $parameterSampleFile2, $parameterSampleFile2arg, $parameterSampleFile2JoinDelimiter ) = get_parameter_sample_files( $config, $section, "parameterSampleFile2" );
  my ( $parameterSampleFile3, $parameterSampleFile3arg, $parameterSampleFile3JoinDelimiter ) = get_parameter_sample_files( $config, $section, "parameterSampleFile3" );

  my ( $scatter, $scatter_arg, $scatter_JoinDelimiter ) = get_parameter_sample_files( $config, $section, "scatter" );
  #print(Dumper($scatter));

  $option = get_parameter_file_option($config, $section, $option);

  my $shfile= $self->get_task_filename( $pbs_dir, $task_name );
  my $sh;
  open( $sh, ">$shfile" ) or die "Cannot create $shfile";
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
      $param_option2 = $parameterSampleFile2arg . " " . basename($file2);
    }
  }

  my $param_option3 = "";
  if ( not $bFound3 ) {
    my $file3 = save_parameter_sample_file( $config, $section, "parameterSampleFile3", "${result_dir}/${task_name}_${task_suffix}_fileList3.list" );
    if ( $file3 ne "" ) {
      $param_option3 = $parameterSampleFile3arg . " " . basename($file3);
    }
  }

  my $expect_result = $self->result( $config, $section );
  #print(Dumper($expect_result));

  for my $sample_name ( sort keys %$parameterSampleFile1 ) {
    my $cur_dir = $output_to_same_folder ? $result_dir : $result_dir . "/$sample_name";
    if(! -e $cur_dir){
      create_directory_or_die($cur_dir);
    }
    for my $scatter_name (sort keys %$scatter){
      my $curname = $sample_name . "_" . $scatter_name;

      my $curOption = $option;

      my $pbs_file = $self->get_pbs_filename( $pbs_dir, $curname );
      my $pbs_name = basename($pbs_file);
      my $log      = $self->get_log_filename( $log_dir, $curname );

      my $log_desc = $cluster->get_log_description($log);

      my $final_file            = $expect_result->{$curname}[-1];
      my $pbs                   = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );


      print $sh "if [[ ! -s $final_file ]]; then
  \$MYCMD ./$pbs_name 
fi

";

      my $final_prefix = $curname . $output_file_prefix;

      my $output_option = "$output_arg $final_prefix";
      if ($curOption =~ /__OUTPUT__/){
        $curOption =~ s/__OUTPUT__/$final_prefix/g;
        $output_option = "";
      }

      if (option_contains_arg($curOption, $output_arg)) {
        $output_option = "";
      }

      if ($curOption =~ /__NAME__/){
        $curOption =~ s/__NAME__/$sample_name/g;
      }

      if ($curOption =~ /__FILE__/){
        my $param_option1 = get_program_param( $parameterSampleFile1, "", $parameterSampleFile1JoinDelimiter, $sample_name, $result_dir, 1);
        $curOption =~ s/__FILE__/$param_option1/g;
      } elsif (option_contains_arg($curOption, $parameterSampleFile1arg)) {
      } else{
        my $param_option1 = get_program_param( $parameterSampleFile1, $parameterSampleFile1arg, $parameterSampleFile1JoinDelimiter, $sample_name, $result_dir, 1 );
        $curOption = $curOption . " " . $param_option1;
      }

      if ($curOption =~ /__SCATTER__/){
        my $param_option1 = get_program_param( $scatter, "", $scatter_JoinDelimiter, $scatter_name, $result_dir, 1);
        $curOption =~ s/__SCATTER__/$param_option1/g;
      } elsif (option_contains_arg($curOption, $scatter_arg)) {
      } else{
        my $param_option1 = get_program_param( $scatter, $scatter_arg, $scatter_JoinDelimiter, $scatter_name, $result_dir, 1 );
        $curOption = $curOption . " " . $param_option1;
      }
      if ($curOption =~ /__SCATTERNAME__/){
        $curOption =~ s/__SCATTERNAME__/$scatter_name/g;
      }

      if ( not $bFound2 ) {
        $curOption = $curOption . " " . $param_option2;
      }else{
        $curOption = $curOption . " " . get_program_param($parameterSampleFile2, $parameterSampleFile2arg, $parameterSampleFile2JoinDelimiter, $sample_name, $result_dir, 2);
      }

      if ( not $bFound3 ) {
        $curOption = $curOption . " " . $param_option3;
      }else{
        $curOption = $curOption . " " . get_program_param($parameterSampleFile3, $parameterSampleFile3arg, $parameterSampleFile3JoinDelimiter, $sample_name, $result_dir, 3);
      }

      print $pbs "
  $interpretor $program $curOption $output_option

  ";
      $self->close_pbs( $pbs, $pbs_file );
    }
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

  my $samplename_in_result = get_option( $config, $section, "samplename_in_result", 1 );

  my ($source_files, $source_file_arg, $source_file_join_delimiter) = get_parameter_sample_files( $config, $section, "source" );
  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder" );
  my $output_exts = get_output_ext_list( $config, $section );

  my ( $scatter, $scatter_arg, $scatter_JoinDelimiter ) = get_parameter_sample_files( $config, $section, "scatter" );

  my $result = {};
  for my $sample_name ( sort keys %$source_files ) {
    my $cur_dir = $output_to_same_folder ? $result_dir : $result_dir . "/$sample_name";
    for my $scatter_name (sort keys %$scatter){
      my $curname = $sample_name . "_" . $scatter_name;

      my $result_files = [];
      foreach my $cur_ext (@$output_exts) {
        my $final_file = $samplename_in_result ? $curname . $cur_ext : $cur_ext;
        push( @$result_files, "${cur_dir}/$final_file" );
      }
      $result->{$curname} = filter_array( $result_files, $pattern );      
    }
  }

  #print(Dumper($result));
  return $result;
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;
  
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  
  my ( $source_files, $source_file_arg, $source_file_join_delimiter ) = get_parameter_sample_files( $config, $section, "source" );
  my ( $scatter, $scatter_arg, $scatter_JoinDelimiter ) = get_parameter_sample_files( $config, $section, "scatter" );

  my $result = {};
  
  for my $sample_name ( sort keys %$source_files ) {
    for my $scatter_name (sort keys %$scatter){
      my $cur_name = $sample_name . "_" . $scatter_name;
      my $pbs_file = $self->get_pbs_filename( $pbs_dir, $cur_name );
      $result->{$cur_name} = $pbs_file;
    }
  }
  
  return ($result);
}

#get current pbs and its dependent sample names
#TODO: consider other ref links
sub get_pbs_source {
  my ( $self, $config, $section ) = @_;
  
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  
  my ( $source_files, $source_file_arg, $source_file_join_delimiter ) = get_parameter_sample_files( $config, $section, "source" );
  my ( $scatter, $scatter_arg, $scatter_JoinDelimiter ) = get_parameter_sample_files( $config, $section, "scatter" );

  my $result = {};
  
  for my $sample_name ( sort keys %$source_files ) {
    for my $scatter_name (sort keys %$scatter){
      my $cur_name = $sample_name . "_" . $scatter_name;
      my $pbs_file = $self->get_pbs_filename( $pbs_dir, $cur_name );
      $result->{$pbs_file} = [$sample_name];
    }
  }
  
  return ($result);
}

#get result sample name and its depedent current pbs
sub get_result_pbs {
  my ( $self, $config, $section ) = @_;
  
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  
  my ( $source_files, $source_file_arg, $source_file_join_delimiter ) = get_parameter_sample_files( $config, $section, "source" );
  my ( $scatter, $scatter_arg, $scatter_JoinDelimiter ) = get_parameter_sample_files( $config, $section, "scatter" );

  my $result = {};

  for my $sample_name ( sort keys %$source_files ) {
    for my $scatter_name (sort keys %$scatter){
      my $cur_name = $sample_name . "_" . $scatter_name;
      my $pbs_file = $self->get_pbs_filename( $pbs_dir, $cur_name );
      $result->{$cur_name} = $pbs_file;
    }
  }
  
  return ($result);
}

1;
