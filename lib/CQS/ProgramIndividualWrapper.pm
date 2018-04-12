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

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  my $task_suffix = get_option( $config, $section, "suffix", "" );
  $self->{_task_suffix} = $task_suffix;

  my $interpretor = get_option( $config, $section, "interpretor", "" );
  my $program = get_option( $config, $section, "program" );
  if ( !File::Spec->file_name_is_absolute($program) ) {
    $program = dirname(__FILE__) . "/$program";
  }
  if ( !( -e $program ) ) {
    die("program $program defined but not exists!");
  }

  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder" );
  my $output_ext            = get_option( $config, $section, "output_ext", "" );
  my $first_file_only       = get_option( $config, $section, "first_file_only", 0 );
  my $output_arg            = get_option( $config, $section, "output_arg", "" );

  my ( $parameterSampleFile1, $parameterSampleFile1arg ) = get_parameter_sample_files( $config, $section, "source" );
  my ( $parameterSampleFile2, $parameterSampleFile2arg ) = get_parameter_sample_files( $config, $section, "parameterSampleFile2" );
  my ( $parameterSampleFile3, $parameterSampleFile3arg ) = get_parameter_sample_files( $config, $section, "parameterSampleFile3" );

  my ( $parameterFile1, $parameterFile1arg ) = get_parameter_file( $config, $section, "parameterFile1" );
  my ( $parameterFile2, $parameterFile2arg ) = get_parameter_file( $config, $section, "parameterFile2" );
  my ( $parameterFile3, $parameterFile3arg ) = get_parameter_file( $config, $section, "parameterFile3" );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %$parameterSampleFile1 ) {
    my $pfiles1 = $parameterSampleFile1->{$sample_name};
    my $idxend = $first_file_only ? 0 : ( scalar(@$pfiles1) - 1 );

    my $cur_dir = $output_to_same_folder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir );

    for my $i ( 0 .. $idxend ) {
      my $pfile1 = $pfiles1->[$i];
      my $final_file = ($first_file_only or (1 == scalar(@$pfiles1))) ? $sample_name . $output_ext : basename($pfile1) . $output_ext;

      my $curOption = "";
      if ( defined $parameterSampleFile2->{$sample_name} ) {
        $curOption = $curOption . " " . $parameterSampleFile2arg . " " . $parameterSampleFile2->{$sample_name}[$i];
      }

      if ( defined $parameterSampleFile3->{$sample_name} ) {
        $curOption = $curOption . " " . $parameterSampleFile3arg . " " . $parameterSampleFile3->{$sample_name}[$i];
      }

      print $pbs "
if [[ ! -s $final_file ]]; then
  $interpretor $program $option $parameterSampleFile1arg $pfile1 $curOption $parameterFile1arg $parameterFile1 $parameterFile2arg $parameterFile2 $parameterFile3arg $parameterFile3 $output_arg $final_file
fi

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  my $task_suffix = get_option( $config, $section, "suffix", "" );
  $self->{_task_suffix} = $task_suffix;

  my $first_file_only = get_option( $config, $section, "first_file_only", 0 );

  my ( $parameterSampleFile1, $parameterSampleFile1arg ) = get_parameter_sample_files( $config, $section, "source" );
  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder" );
  my $output_ext            = get_option( $config, $section, "output_ext", "" );
  my $output_other_ext      = get_option( $config, $section, "output_other_ext", "" );
  my @output_other_exts;
  if ( $output_other_ext ne "" ) {
    @output_other_exts = split( ",", $output_other_ext );
  }

  my $result = {};
  for my $sample_name ( sort keys %$parameterSampleFile1 ) {
    my $pfiles1 = $parameterSampleFile1->{$sample_name};
    my $idxend = $first_file_only ? 0 : ( scalar(@$pfiles1) - 1 );

    my $cur_dir = $output_to_same_folder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );

    my @result_files = ();
    for my $i ( 0 .. $idxend ) {
      my $pfile1 = $pfiles1->[$i];
      my $final_file = ($first_file_only or (1 == scalar(@$pfiles1))) ? $sample_name . $output_ext : basename($pfile1) . $output_ext;

      push( @result_files, "${cur_dir}/$final_file" );
      if ( $output_other_ext ne "" ) {
        foreach my $output_other_ext_each (@output_other_exts) {
          my $other_file = $first_file_only ? $sample_name . $output_other_ext_each : basename($pfile1) . $output_other_ext_each;
          push( @result_files, "${cur_dir}/$other_file" );
        }
      }
    }

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
