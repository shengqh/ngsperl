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
use GATK4::GATK4UniqueTask;

our @ISA = qw(GATK4::GATK4UniqueTask);

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $java_option = $self->get_java_option( $config, $section, $memory );

  #parameter files
  $self->get_docker_value(1);

  my $intervals = parse_param_file( $config, $section, "preprocessed_intervals", 1 );
  my $blacklist_intervals = get_param_file( $config->{$section}{blacklist_file}, "blacklist_file", 0 );
  my $blacklist_intervals_option = $blacklist_intervals ? "-XL " . $blacklist_intervals : "";

  my $ref_fasta_dict = get_param_file( $config->{$section}{ref_fasta_dict}, "ref_fasta_dict", 1 );
  my $ref_fasta      = get_param_file( $config->{$section}{ref_fasta},      "ref_fasta",      1 );

  my $contig_ploidy_priors = get_param_file($config->{$section}{contig_ploidy_priors}, "contig_ploidy_priors_file", 0 );
  my $contig_ploidy_priors_option = $contig_ploidy_priors? "-c " . $contig_ploidy_priors:"";

  my $parameters = get_parameter_options(
    $config, $section, "--",
    [
      "minimum-gc-content",                      "maximum-gc-content",                         #
      "minimum-mappability",                     "maximum-mappability",                        #
      "minimum-segmental-duplication-content",   "maximum-segmental-duplication-content",      #
      "low-count-filter-count-threshold",        "low-count-filter-percentage-of-samples",     #
      "extreme-count-filter-minimum-percentile", "extreme-count-filter-maximum-percentile",    #
      "extreme-count-filter-percentage-of-samples"
    ],
    [
      "0.1", "0.9",                                                                            #
      "0.9", "1.0",                                                                            #
      "0.0", "0.5",                                                                            #
      "5",   "90.0",                                                                           #
      "1",   "99.0",                                                                           #
      "90.0"
    ]
  );

  my $script = dirname(__FILE__) . "/filterIntervals.py";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $gc_file = $task_name . ".annotated.tsv";
  my $final_file = $task_name . ".filtered.interval_list";
  my $final_file_tmp = $task_name . ".filtered.all.interval_list";
  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $inputOption = "";
  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sampleFile   = $sample_files[0];
    $inputOption = $inputOption . " \\\n  --input " . $sampleFile;
  }

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command );
  print $pbs "  
export HOME=$result_dir
export PYTHONPATH=

source activate gatk

cd $result_dir
gatk --java-options \"$java_option\" AnnotateIntervals \\
  -L ${intervals} $blacklist_intervals_option \\
  -R $ref_fasta \\
  -imr OVERLAPPING_ONLY \\
  -O $gc_file

gatk --java-options \"$java_option\" FilterIntervals $option \\
  -L ${intervals} $blacklist_intervals_option $inputOption \\
  --annotated-intervals $gc_file \\
  --interval-merging-rule OVERLAPPING_ONLY $parameters \\
  --output $final_file_tmp

python3 $script -i $final_file_tmp -o $final_file $contig_ploidy_priors_option

#rm $final_file_tmp
";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $final_file = $task_name . ".filtered.interval_list";
  my $result = { $task_name => filter_array( ["${result_dir}/${final_file}"], $pattern ) };

  return $result;
}

1;
