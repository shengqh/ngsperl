#!/usr/bin/perl
package GATK4::PreprocessIntervals;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::UniqueTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_pi";
  bless $self, $class;
  return $self;
}

#Based on https://github.com/broadinstitute/gatk/blob/master/scripts/cnv_wdl/cnv_common_tasks.wdl
sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = get_parameter( $config, $section );

  my $java_option = $self->get_java_option( $config, $section, $memory );

  #parameter files
  #my $gatkJar = get_param_file( $config->{$section}{gatk_jar}, "gatk_jar", 1 );
  my $gatk4_singularity = get_param_file( $config->{$section}{gatk4_singularity}, "gatk4_singularity", 1 );

  my $intervals           = get_param_file( $config->{$section}{interval_file},  "interval_file",  1 );
  my $blacklist_intervals = get_param_file( $config->{$section}{blacklist_file}, "blacklist_file", 0 );
  my $blacklist_intervals_option = $blacklist_intervals ? "-XL " . $blacklist_intervals : "";

  my $ref_fasta_dict = get_param_file( $config->{$section}{ref_fasta_dict}, "ref_fasta_dict", 1 );
  my $ref_fasta      = get_param_file( $config->{$section}{ref_fasta},      "ref_fasta",      1 );

  #use default value in software rather than assign here since GATK team is still tunning the parameters.
  my $parameters = get_parameter_options( $config, $section, "--", [ "padding", "bin-length" ] );

  my $final_file = $task_name . ".preprocessed.interval_list";

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh "  
export HOME=$result_dir
export PYTHONPATH=

source activate gatk

cd $result_dir

gatk --java-options \"$java_option\" PreprocessIntervals \\
  -L $intervals $blacklist_intervals_option \\
  --sequence-dictionary $ref_fasta_dict \\
  --reference $ref_fasta \\
  --interval-merging-rule OVERLAPPING_ONLY $parameters \\
  --output $final_file
  
rm -rf .conda

";
  close($sh);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command );
  print $pbs "singularity exec $gatk4_singularity bash $shfile \n";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $final_file = $task_name . ".preprocessed.interval_list";
  my $result = { $task_name => filter_array( ["${result_dir}/${final_file}"], $pattern ) };
  return $result;
}

1;
