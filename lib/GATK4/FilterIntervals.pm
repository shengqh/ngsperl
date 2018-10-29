#!/usr/bin/perl
package GATK4::FilterIntervals;

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
  $self->{_suffix} = "_fi";
  bless $self, $class;
  return $self;
}

#Based on https://github.com/broadinstitute/gatk/blob/master/scripts/cnv_wdl/cnv_common_tasks.wdl
sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = get_parameter( $config, $section );

  my $java_option = $self->get_java_option($config, $section, $memory);

  #parameter files
  my $gatk_singularity = get_param_file( $config->{$section}{gatk_singularity}, "gatk_singularity", 1 );

  my $intervals = parse_param_file( $config, $section, "preprocessed_intervals", 1 );
  my $blacklist_intervals = get_param_file( $config->{$section}{blacklist_file}, "blacklist_file", 0 );
  my $blacklist_intervals_option = $blacklist_intervals ? "-XL " . $blacklist_intervals : "";

  my $ref_fasta_dict = get_param_file( $config->{$section}{ref_fasta_dict}, "ref_fasta_dict", 1 );
  my $ref_fasta      = get_param_file( $config->{$section}{ref_fasta},      "ref_fasta",      1 );

  my $minimum_gc_content                         = get_option( $config, $section, "minimum_gc_content",                         0.1 );
  my $maximum_gc_content                         = get_option( $config, $section, "maximum_gc_content",                         0.9 );
  my $minimum_mappability                        = get_option( $config, $section, "minimum_mappability",                        0.9 );
  my $maximum_mappability                        = get_option( $config, $section, "maximum_mappability",                        1.0 );
  my $minimum_segmental_duplication_content      = get_option( $config, $section, "minimum_segmental_duplication_content",      0.0 );
  my $maximum_segmental_duplication_content      = get_option( $config, $section, "maximum_segmental_duplication_content",      0.5 );
  my $low_count_filter_count_threshold           = get_option( $config, $section, "low_count_filter_count_threshold",           5 );
  my $low_count_filter_percentage_of_samples     = get_option( $config, $section, "low_count_filter_percentage_of_samples",     90.0 );
  my $extreme_count_filter_minimum_percentile    = get_option( $config, $section, "extreme_count_filter_minimum_percentile",    1.0 );
  my $extreme_count_filter_maximum_percentile    = get_option( $config, $section, "extreme_count_filter_maximum_percentile",    99.0 );
  my $extreme_count_filter_percentage_of_samples = get_option( $config, $section, "extreme_count_filter_percentage_of_samples", 90.0 );

  my $final_file = $task_name . ".filtered.interval_list";
  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $inputOption = "";
  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sampleFile   = $sample_files[0];
    $inputOption = $inputOption . " --input " . $sampleFile;
  }


  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh "  
gatk --java-options \"$java_option\" FilterIntervals \\
            -L ${intervals} $blacklist_intervals_option \\
            $inputOption \\
            --minimum-gc-content ${minimum_gc_content} \\
            --maximum-gc-content ${maximum_gc_content} \\
            --minimum-mappability ${minimum_mappability} \\
            --maximum-mappability ${maximum_mappability} \\
            --minimum-segmental-duplication-content ${minimum_segmental_duplication_content} \\
            --maximum-segmental-duplication-content ${maximum_segmental_duplication_content} \\
            --low-count-filter-count-threshold ${low_count_filter_count_threshold} \\
            --low-count-filter-percentage-of-samples ${low_count_filter_percentage_of_samples} \\
            --extreme-count-filter-minimum-percentile ${extreme_count_filter_minimum_percentile} \\
            --extreme-count-filter-maximum-percentile ${extreme_count_filter_maximum_percentile} \\
            --extreme-count-filter-percentage-of-samples ${extreme_count_filter_percentage_of_samples} \\
            --interval-merging-rule OVERLAPPING_ONLY \\
            --output $final_file
";
  close($sh);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command );
  print $pbs "singularity run $gatk_singularity $shfile \n";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $final_file = $task_name . ".filtered.interval_list";
  my $result = { $task_name => filter_array( ["${result_dir}/${final_file}"], $pattern ) };

  return $result;
}

1;
