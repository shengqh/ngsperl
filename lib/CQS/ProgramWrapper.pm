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
use Data::Dumper;

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
  my $input;
  ($cur_option, $input) = process_parameter_sample_file($config, $section, $result_dir, $task_name, $task_suffix, $cur_option, $source_key, 1);

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

  my $cur_post_command = get_option( $config, $section, "post_command", "" );
  if ($cur_post_command =~ /__NAME__/){
    $cur_post_command =~ s/__NAME__/$sample_name/g;
  }

  if ($cur_post_command =~ /__FILE__/){
    $cur_post_command =~ s/__FILE__/$input/g;
  }

  if ($cur_post_command =~ /__OUTPUT__/){
    $cur_post_command =~ s/__OUTPUT__/$final_prefix/g;
  }

  return($cur_option, $output_option, $cur_init_command, $cur_post_command);
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  my $task_suffix = get_option( $config, $section, "suffix", "" );
  $self->{_task_suffix} = $task_suffix;

  my $check_file_ext = get_option($config, $section, "check_file_ext", "");

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

  my $source_key = "source";
  if ((not defined $config->{$section}{"source"} ) && (not defined $config->{$section}{"source_ref"} )) {
    $source_key = "parameterSampleFile1";
  }

  for my $index (2..10){
    my $key = "parameterSampleFile" . $index;
    if (not has_raw_files($config, $section, $key)){
      next;
    }

    my ( $parameterSampleFile, $parameterSampleFilearg, $parameterSampleFileJoinDelimiter ) = get_parameter_sample_files( $config, $section, $key );
    my $listfile = save_parameter_sample_file( $config, $section, $key, "${result_dir}/${task_name}_${task_suffix}_fileList${index}.list" );
    if ( $listfile ne "" ) {
      $listfile = basename($listfile);
      $option = $option . " " . $parameterSampleFilearg . " " . $listfile;
    }
  }

  $option = $option . " " . get_parameter_file_option($config, $section);

  my $pbs_file   = $self->get_pbs_filename( $pbs_dir, $task_name, ".pbs" );
  my $pbs_name   = basename($pbs_file);
  my $log        = $self->get_log_filename( $log_dir, $task_name, ".log" );
  my $log_desc   = $cluster->get_log_description($log);

  my $results = $self->result( $config, $section );
  my $result_files = $results->{$task_name};

  #print(Dumper($result_files));
  
  my $pbs;
  my $final_file;
  if (defined $result_files){
    $final_file = $result_files->[0];

    if ($final_file ne $result_dir) {
      if (rindex($final_file, $result_dir, 0) == 0){
        $final_file=substr($final_file, length($result_dir)+1);
      }
    }
    
    my $checkFile = $result_files->[-1] . $check_file_ext;
    $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $checkFile );
  }else{
    $final_file = $target_dir . "/";
    if(defined $results){
      my @sampleKeys = reverse sort(keys (%$results));
      my $firstKey = $sampleKeys[0];
      my $checkFile = $results->{$firstKey}[-1] . $check_file_ext;
      $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $checkFile);
    }else{
      $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir);
    }
  }

  #print($final_file);

  my ($cur_option, $output_option, $cur_init_command, $cur_post_command) = replace_tag( $config, $section, $task_name, $result_dir, $option, $source_key, $task_name, $final_file);

  print $pbs "
$cur_init_command

$interpretor $program $cur_option $output_option

$cur_post_command
";

  $self->close_pbs( $pbs, $pbs_file );

  print "!!!pbs file $pbs_file created, you can run this pbs file to submit to cluster.\n";
}

1;
