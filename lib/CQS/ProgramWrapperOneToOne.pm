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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  my $task_suffix = get_option( $config, $section, "suffix", "" );
  $self->{_task_suffix} = $task_suffix;

  if ($option =~ /__MEMORY__/) {
    $option =~ s/__MEMORY__/$memory/g;
  }

  my $init_command = get_option( $config, $section, "init_command", "" );

  my $interpretor = get_option( $config, $section, "interpretor", "" );
  my $program     = get_program( $config, $section );

  my $output_to_folder = get_option( $config, $section, "output_to_folder", 0 );
  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder", 1);
  my $output_file_prefix    = get_option( $config, $section, "output_file_prefix", (!$output_to_folder) );
  my $output_arg            = get_option( $config, $section, "output_arg" );
  my $no_output            = get_option( $config, $section, "no_output", 0 );
  
  my ( $parameterSampleFile1, $parameterSampleFile1arg, $parameterSampleFile1JoinDelimiter ) = get_parameter_sample_files( $config, $section, "source" );
  my @sample_names = ( sort keys %$parameterSampleFile1 );
  my $has_multi_samples = scalar(@sample_names) > 1;

  $option = $option . " " . get_parameter_file_option($config, $section);

  my $shfile;
  my $sh;
  if ($has_multi_samples) {
    $shfile = $self->get_task_filename( $pbs_dir, $task_name );
    open( $sh, ">$shfile" ) or die "Cannot create $shfile";
    print $sh get_run_command($sh_direct) . "\n";
  }

  my $paramFileMap = {};
  for my $index (2..10){
    my $key = "parameterSampleFile" . $index;
    if (not has_raw_files($config, $section, $key)){
      next;
    }

    my ( $parameterSampleFile, $parameterSampleFilearg, $parameterSampleFileJoinDelimiter, $paramterSampleFileNameArg, $parameterSampleFileNameJoinDelimiter ) = get_parameter_sample_files( $config, $section, $key );
    my $bFound = 0;
    for my $sample_name ( sort keys %$parameterSampleFile1 ) {
      if ( defined $parameterSampleFile->{$sample_name} ) {
        $bFound = 1;
      }
    }
    if ( not $bFound ) {
      my $input;
      ($option, $input) = process_parameter_sample_file($config, $section, $result_dir, $task_name, $task_suffix, $option, $key, $index);
    }else{
      $paramFileMap->{$key} = [$parameterSampleFile, $parameterSampleFilearg, $parameterSampleFileJoinDelimiter, $index];
    }
  }

  my $expect_result = $self->result( $config, $section );

  for my $sample_name ( sort keys %$parameterSampleFile1 ) {
    my $curOption = $option;

    my $cur_dir = $output_to_same_folder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $log_desc = $cluster->get_log_description($log);

    my $final_file = $expect_result->{$sample_name}[-1];
    my $pbs        = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    if ( $has_multi_samples ) {
      print $sh "if [[ ! -s $final_file ]]; then
  \$MYCMD ./$pbs_name 
fi
";
    }

    my $final_prefix = $output_to_folder ? "." : $sample_name . $output_file_prefix;

    my $localized_files = [];
    $parameterSampleFile1->{$sample_name} = $self->localize_files_in_tmp_folder($pbs, $parameterSampleFile1->{$sample_name}, $localized_files);

    #print Dumper($parameterSampleFile1->{$sample_name});

    if ($curOption =~ /__NAME__/){
      $curOption =~ s/__NAME__/$sample_name/g;
    }

    if ($curOption =~ /__FILE__/){
      my $param_option1 = get_program_param( $parameterSampleFile1, "", $parameterSampleFile1JoinDelimiter, $sample_name, $result_dir, 1 );
      #print("delimiter=" . $parameterSampleFile1JoinDelimiter . "\n");
      $curOption =~ s/__FILE__/$param_option1/g;
    } elsif (option_contains_arg($curOption, $parameterSampleFile1arg)) {
    } else{
      my $param_option1 = get_program_param( $parameterSampleFile1, $parameterSampleFile1arg, $parameterSampleFile1JoinDelimiter, $sample_name, $result_dir, 1 );
      $curOption = $curOption . " " . $param_option1;
    }

    my $ignored_map = {};
    for my $index (2..10){
      my $key = "parameterSampleFile" . $index;
      $ignored_map->{$key} = 0;
      my $place_hold = "__FILE${index}__";
      if ($curOption =~ /$place_hold/){
        my $file_options = $paramFileMap->{$key};
        #$parameterSampleFile, $parameterSampleFilearg, $parameterSampleFileJoinDelimiter, $index
        my $param_option = get_program_param( $file_options->[0], "", $file_options->[2], $sample_name, $result_dir, $index );
        #print("delimiter=" . $parameterSampleFile1JoinDelimiter . "\n");
        $curOption =~ s/$place_hold/$param_option/g;
        $ignored_map->{$key} = 1;
      }
    }

    my $output_option = "$output_arg $final_prefix";
    if ($curOption =~ /__OUTPUT__/){
      $curOption =~ s/__OUTPUT__/$final_prefix/g;
      $output_option = "";
    }

    if (option_contains_arg($curOption, $output_arg)) {
      $output_option = "";
    }

    if ($no_output){
      $output_option = "";
    }

    my $cur_init_command = $init_command;
    if ($cur_init_command =~ /__NAME__/){
      $cur_init_command =~ s/__NAME__/$sample_name/g;
    }

    if ($cur_init_command =~ /__FILE__/){
      my $param_option1 = get_program_param( $parameterSampleFile1, "", $parameterSampleFile1JoinDelimiter, $sample_name, $result_dir, 1 );
      $cur_init_command =~ s/__FILE__/$param_option1/g;
    }

    for my $key (sort keys %$paramFileMap){
      my $values = $paramFileMap->{$key};
      if(!$ignored_map->{$key}) {
        $curOption = $curOption . " " . get_program_param( $values->[0], $values->[1],$values->[2], $sample_name, $result_dir, $values->[3] );
      }
    }

    print $pbs "
$cur_init_command    

$interpretor $program $curOption $output_option

";

    $self->clean_temp_files($pbs, $localized_files);

    $self->close_pbs( $pbs, $pbs_file );
  }

  if ( $has_multi_samples ) {
    close $sh;
    if ( is_linux() ) {
      chmod 0755, $shfile;
    }

    print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my ( $source_files, $source_file_arg, $source_file_join_delimiter ) = get_parameter_sample_files( $config, $section, "source" );
  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder" );
  my $output_exts           = get_output_ext_list( $config, $section );

  my $result = {};
  for my $sample_name ( sort keys %$source_files ) {
    my $cur_dir = $output_to_same_folder ? $result_dir : $result_dir . "/$sample_name";

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
