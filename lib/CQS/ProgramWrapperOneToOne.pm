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
  $self->{_name}                  = __PACKAGE__;
  $self->{_suffix}                = "_o2o";
  $self->{_output_to_same_folder} = 1;
  bless $self, $class;
  return $self;
} ## end sub new


sub get_pbs_key {
  my ( $self, $config, $section ) = @_;

  my $has_source               = has_raw_files( $config, $section, "source" );
  my $has_parameterSampleFile1 = has_raw_files( $config, $section, "parameterSampleFile1" );

  if ( $has_source && $has_parameterSampleFile1 ) {
    die("Cannot define source and parameterSampleFile1 at same time for section $section");
  }

  if ($has_source) {
    return ("source");
  }

  if ($has_parameterSampleFile1) {
    return ("parameterSampleFile1");
  }

  die("Define source or parameterSampleFile1 for section $section first!");
} ## end sub get_pbs_key


sub print_sh_pbs {
  my ( $self, $sh, $pbs_name, $final_file, $pbs_index, $cur_sh_log ) = @_;
  print $sh "if [[ ! -s $final_file ]]; then
  \$MYCMD ./$pbs_name $cur_sh_log
fi
";
} ## end sub print_sh_pbs


sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command, $sh_log ) = $self->init_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  my $task_suffix = get_option( $config, $section, "suffix", "" );
  $self->{_task_suffix} = $task_suffix;

  if ( $option =~ /__MEMORY__/ ) {
    $option =~ s/__MEMORY__/$memory/g;
  }

  my $check_file_ext = get_option( $config, $section, "check_file_ext", "" );

  my $post_command = get_option( $config, $section, "post_command", "" );

  my $interpretor = get_option( $config, $section, "interpretor", "" );
  my $program     = get_program( $config, $section );

  my $output_to_folder      = get_option( $config, $section, "output_to_folder",      0 );
  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder", $self->{_output_to_same_folder} );
  my $output_file_prefix    = get_option( $config, $section, "output_file_prefix",    "" );
  my $output_arg            = get_option( $config, $section, "output_arg",            "" );
  my $no_output             = get_option( $config, $section, "no_output",             0 );
  my $no_input              = get_option( $config, $section, "no_input",              0 );

  my $list_files = get_option( $config, $section, "copy_files", "" );

  my $other_localization_ext_array = get_option( $config, $section, "other_localization_ext_array", [] );

  my $source_key = $self->get_pbs_key( $config, $section );
  my ( $parameterSampleFile1, $parameterSampleFile1arg, $parameterSampleFile1JoinDelimiter ) = get_parameter_sample_files( $config, $section, $source_key );

  my @sample_names      = ( sort keys %$parameterSampleFile1 );
  my $has_multi_samples = scalar(@sample_names) > 1;

  $option = get_parameter_file_option( $config, $section, $option );

  my $shfile;
  my $sh;
  if ($has_multi_samples) {
    $shfile = $self->get_task_filename( $pbs_dir, $task_name );
    open( $sh, ">$shfile" ) or die "Cannot create $shfile";
    print $sh get_run_command($sh_direct) . "\n";
  }

  my $paramFileMap = {};
  for my $index ( 2 .. 10 ) {
    my $key = "parameterSampleFile" . $index;
    if ( not has_raw_files( $config, $section, $key ) ) {
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
      ( $option, $input ) = process_parameter_sample_file( $config, $section, $result_dir, $task_name, $task_suffix, $option, $key, $index );
    }
    else {
      $paramFileMap->{$key} = [ $parameterSampleFile, $parameterSampleFilearg, $parameterSampleFileJoinDelimiter, $index ];
    }
  } ## end for my $index ( 2 .. 10)

  my $expect_result = $self->get_expect_result_for_perform( $config, $section );

  my $can_result_be_empty_file = get_option( $config, $section, "can_result_be_empty_file", 0 );

  my $pbs_index = 0;
  for my $sample_name ( sort keys %$parameterSampleFile1 ) {
    $pbs_index = $pbs_index + 1;
    my $curOption = $option;

    my $cur_dir = $output_to_same_folder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );
    copy_files( $list_files, $cur_dir );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $log_desc = $cluster->get_log_description($log);

    my $final_file = $check_file_ext ne "" ? $sample_name . $check_file_ext : $expect_result->{$sample_name}[-1];
    #print("final file=" . $final_file . "\n");
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file, $init_command, $can_result_be_empty_file );

    my $cur_sh_log = $sh_log;
    if ( $cur_sh_log =~ /__NAME__/ ) {
      $cur_sh_log =~ s/__NAME__/$sample_name/g;
    }

    if ($has_multi_samples) {
      $self->print_sh_pbs( $sh, $pbs_name, $final_file, $pbs_index, $cur_sh_log );
    }

    my $final_prefix = $output_to_folder ? "." : $sample_name . $output_file_prefix;

    my $localized_files = [];
    $parameterSampleFile1->{$sample_name} = $self->localize_files_in_tmp_folder( $pbs, $parameterSampleFile1->{$sample_name}, $localized_files, $other_localization_ext_array, get_option( $config, $section, "no_bai", 0 ) );

    #print Dumper($parameterSampleFile1->{$sample_name});

    if ( $curOption =~ /__NAME__/ ) {
      $curOption =~ s/__NAME__/$sample_name/g;
    }

    if ( !$no_input ) {
      if ( $curOption =~ /__FILE__/ ) {
        my $param_option1 = get_program_param( $parameterSampleFile1, "", $parameterSampleFile1JoinDelimiter, $sample_name, $result_dir, 1 );
        #print("delimiter=" . $parameterSampleFile1JoinDelimiter . "\n");
        $curOption =~ s/__FILE__/$param_option1/g;
      }
      elsif ( option_contains_arg( $curOption, $parameterSampleFile1arg ) ) {
      }
      else {
        my $param_option1 = get_program_param( $parameterSampleFile1, $parameterSampleFile1arg, $parameterSampleFile1JoinDelimiter, $sample_name, $result_dir, 1 );
        $curOption = $curOption . " " . $param_option1;
      }
    } ## end if ( !$no_input )

    my $ignored_map = {};
    for my $index ( 2 .. 10 ) {
      my $key = "parameterSampleFile" . $index;
      $ignored_map->{$key} = 0;
      my $place_hold = "__FILE${index}__";
      if ( $curOption =~ /$place_hold/ ) {
        my $file_options = $paramFileMap->{$key};
        #$parameterSampleFile, $parameterSampleFilearg, $parameterSampleFileJoinDelimiter, $index
        my $param_option = get_program_param( $file_options->[0], "", $file_options->[2], $sample_name, $result_dir, $index );
        #print("delimiter=" . $parameterSampleFile1JoinDelimiter . "\n");
        $curOption =~ s/$place_hold/$param_option/g;
        $ignored_map->{$key} = 1;
      } ## end if ( $curOption =~ /$place_hold/)
    } ## end for my $index ( 2 .. 10)

    my $output_option = "$output_arg $final_prefix";
    if ( $curOption =~ /__OUTPUT__/ ) {
      $curOption =~ s/__OUTPUT__/$final_prefix/g;
      $output_option = "";
    }

    if ( option_contains_arg( $curOption, $output_arg ) ) {
      $output_option = "";
    }

    if ($no_output) {
      $output_option = "";
    }

    #print("output_option=" . $output_option. "\n");

    my $cur_init_command = $init_command;
    if ( $cur_init_command =~ /__NAME__/ ) {
      $cur_init_command =~ s/__NAME__/$sample_name/g;
    }

    if ( $cur_init_command =~ /__FILE__/ ) {
      my $param_option1 = get_program_param( $parameterSampleFile1, "", $parameterSampleFile1JoinDelimiter, $sample_name, $result_dir, 1 );
      $cur_init_command =~ s/__FILE__/$param_option1/g;
    }

    my $cur_post_command = $post_command;
    if ( $cur_post_command =~ /__NAME__/ ) {
      $cur_post_command =~ s/__NAME__/$sample_name/g;
    }

    if ( $cur_post_command =~ /__FILE__/ ) {
      my $param_option1 = get_program_param( $parameterSampleFile1, "", $parameterSampleFile1JoinDelimiter, $sample_name, $result_dir, 1 );
      $cur_post_command =~ s/__FILE__/$param_option1/g;
    }

    for my $key ( sort keys %$paramFileMap ) {
      my $values = $paramFileMap->{$key};
      if ( !$ignored_map->{$key} ) {
        $curOption = $curOption . " " . get_program_param( $values->[0], $values->[1], $values->[2], $sample_name, $result_dir, $values->[3] );
      }
    } ## end for my $key ( sort keys...)

    print $pbs "
$cur_init_command    

$interpretor $program $curOption $output_option

$cur_post_command
";

    $self->clean_temp_files( $pbs, $localized_files );

    $self->close_pbs( $pbs, $pbs_file );
  } ## end for my $sample_name ( sort...)

  if ($has_multi_samples) {
    close $sh;
    if ( is_linux() ) {
      chmod 0755, $shfile;
    }

    print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
  } ## end if ($has_multi_samples)
} ## end sub perform


sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $samplename_in_result = get_option( $config, $section, "samplename_in_result", 1 );

  my $source_key = $self->get_pbs_key( $config, $section );
  my ( $source_files, $source_file_arg, $source_file_join_delimiter ) = get_parameter_sample_files( $config, $section, $source_key );

  my $output_to_same_folder         = get_option( $config, $section, "output_to_same_folder", $self->{_output_to_same_folder} );
  my $output_exts                   = get_output_ext_list( $config, $section );
  my $output_by_file                = get_option( $config, $section, "output_by_file",                0 );
  my $output_by_file_remove_pattern = get_option( $config, $section, "output_by_file_remove_pattern", "" );

  my $result = {};
  for my $sample_name ( sort keys %$source_files ) {
    my $cur_dir   = $output_to_same_folder ? $result_dir : $result_dir . "/$sample_name";
    my $cur_files = $source_files->{$sample_name};

    my @result_files = ();
    for my $output_ext (@$output_exts) {
      my $result_file;
      if ( $output_ext ne "" ) {
        if ($output_by_file) {
          for my $cur_file (@$cur_files) {
            my $cur_name = basename($cur_file);
            if ( $output_by_file_remove_pattern ne "" ) {
              $cur_name =~ s/$output_by_file_remove_pattern//g;
            }
            push( @result_files, "${cur_dir}/${cur_name}${output_ext}" );
          } ## end for my $cur_file (@$cur_files)
          next;
        } ## end if ($output_by_file)

        if ( $output_ext =~ /__NAME__/ ) {
          $result_file = $output_ext;
          $result_file =~ s/__NAME__/$sample_name/g;
        }
        elsif ( not $samplename_in_result ) {
          $result_file = $output_ext;
        }
        else {
          $result_file = $sample_name . $output_ext;
        }
        push( @result_files, "${cur_dir}/$result_file" );
      } ## end if ( $output_ext ne "")
    } ## end for my $output_ext (@$output_exts)

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  } ## end for my $sample_name ( sort...)
  return $result;
} ## end sub result

1;
