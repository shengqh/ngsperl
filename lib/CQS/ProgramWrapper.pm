#!/usr/bin/perl
package CQS::ProgramWrapper;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::UniqueWrapper;
use File::Spec;

our @ISA = qw(CQS::UniqueWrapper);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_pw";
  bless $self, $class;
  return $self;
}

sub replace_tag {
  my ($config, $section, $task_name, $result_dir, $option, $source_key, $sample_name, $final_prefix) = @_;

  my $task_suffix = get_option( $config, $section, "suffix", "" );

  my $cur_option = $option;
  if ($cur_option =~ /__NAME__/){
    $cur_option =~ s/__NAME__/$sample_name/g;
  }

  my ( $parameterSampleFile1, $parameterSampleFile1arg, $parameterSampleFile1JoinDelimiter ) = get_parameter_sample_files( $config, $section, $source_key );
  my $input = "";
  if (defined $config->{$section}{$source_key . "_type"} && ($config->{$section}{$source_key . "_type"} eq "array")){
    my ( $parameterSampleFile1, $parameterSampleFile1arg, $parameterSampleFile1JoinDelimiter ) = get_parameter_sample_files( $config, $section, "source" );
    $input = get_joined_files($parameterSampleFile1, $parameterSampleFile1JoinDelimiter);
  }else{
    my $parameterSampleFile1 = save_parameter_sample_file( $config, $section, $source_key, "${result_dir}/${task_name}_${task_suffix}_fileList1.list" );
    if($parameterSampleFile1 ne ""){
      $input = basename($parameterSampleFile1);
    }
  }

  if ($cur_option =~ /__FILE__/){
    $cur_option =~ s/__FILE__/$input/g;
  } elsif (option_contains_arg($cur_option, $parameterSampleFile1arg)) {
  } else{
    my $param_option1 = get_program_param( $parameterSampleFile1, $parameterSampleFile1arg, $parameterSampleFile1JoinDelimiter, $sample_name );
    $cur_option = $cur_option . " " . $parameterSampleFile1arg . " " . $input;
  }

  my $output_arg            = get_option( $config, $section, "output_arg" );
  my $no_output            = get_option( $config, $section, "no_output", 0 );

  my $output_option = "$output_arg $final_prefix";
  if ($cur_option =~ /__OUTPUT__/){
    $cur_option =~ s/__OUTPUT__/$final_prefix/g;
    $output_option = "";
  }

  if (option_contains_arg($cur_option, $output_arg)) {
    $output_option = "";
  }

  if ($no_output){
    $output_option = "";
  }

  my $cur_init_command = get_option( $config, $section, "init_command", "" );
  if ($cur_init_command =~ /__NAME__/){
    $cur_init_command =~ s/__NAME__/$sample_name/g;
  }

  if ($cur_init_command =~ /__FILE__/){
    $cur_init_command =~ s/__FILE__/$input/g;
  }

  return($cur_option, $output_option, $cur_init_command);
}

sub get_joined_files {
  my ( $parameterSampleFile1, $join_delimiter ) = @_;
  my $pfiles                  = [];
  for my $individual_sample_name (sort keys %$parameterSampleFile1) {
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
  my $program = get_option( $config, $section, "program" );

  if (get_option($config, $section, "check_program", 1)){
    if ( !File::Spec->file_name_is_absolute($program) ) {
      $program = dirname(__FILE__) . "/$program";
    }
    if ( !( -e $program ) ) {
      die("program $program defined but not exists!");
    }
  }

  my $output_ext = get_option( $config, $section, "output_ext", "" );
  if($output_ext eq ""){
     $output_ext = get_option( $config, $section, "output_file_ext", "" );
  }

  my $source_key = "source";
  if ((not defined $config->{$section}{"source"} ) && (not defined $config->{$section}{"source_ref"} )) {
    $source_key = "parameterSampleFile1";
  }

  my $parameterSampleFile2 = save_parameter_sample_file( $config, $section, "parameterSampleFile2", "${result_dir}/${task_name}_${task_suffix}_fileList2.list" );
  if($parameterSampleFile2 ne ""){
    $parameterSampleFile2 = basename($parameterSampleFile2);
  }
  my $parameterSampleFile2arg = get_option($config, $section, "parameterSampleFile2_arg", "");

  my $parameterSampleFile3 = save_parameter_sample_file( $config, $section, "parameterSampleFile3", "${result_dir}/${task_name}_${task_suffix}_fileList3.list" );
  if($parameterSampleFile3 ne ""){
    $parameterSampleFile3 = basename($parameterSampleFile3);
  }
  my $parameterSampleFile3arg = get_option($config, $section, "parameterSampleFile3_arg", "");

  my $parameterFile1 = parse_param_file( $config, $section, "parameterFile1", 0 );
  my $parameterFile1arg = get_option($config, $section, "parameterFile1_arg", "");

  my $parameterFile2 = parse_param_file( $config, $section, "parameterFile2", 0 );
  my $parameterFile2arg = get_option($config, $section, "parameterFile2_arg", "");

  my $parameterFile3 = parse_param_file( $config, $section, "parameterFile3", 0 );
  my $parameterFile3arg = get_option($config, $section, "parameterFile3_arg", "");

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
  my $log_desc   = $cluster->get_log_description($log);

  my $results = $self->result( $config, $section );
  my $result_files = $results->{$task_name};
  
  my $pbs;
  my $final_file;
  if (defined $result_files){
    $final_file = $result_files->[-1];
    $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
  }else{
    $final_file = $target_dir . "/";
    
    my @sampleKeys = reverse sort(keys (%$results));
    my $firstKey = $sampleKeys[0];
    my $checkFile = $results->{$firstKey}[-1];
    $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $checkFile);
  }

  my ($cur_option, $output_option, $cur_init_command) = replace_tag( $config, $section, $task_name, $result_dir, $option, $source_key, $task_name, $final_file);

  print $pbs "
$cur_init_command

$interpretor $program $cur_option $parameterSampleFile2arg $parameterSampleFile2 $parameterSampleFile3arg $parameterSampleFile3 $parameterFile1arg $parameterFile1 $parameterFile2arg $parameterFile2 $parameterFile3arg $parameterFile3 $output_option
";

  $self->close_pbs( $pbs, $pbs_file );

  print "!!!pbs file $pbs_file created, you can run this pbs file to submit to cluster.\n";
}

1;
