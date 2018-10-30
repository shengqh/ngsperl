#!/usr/bin/perl
package GATK4::GermlineCNVCaller;

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
  $self->{_suffix} = "_gcc";
  bless $self, $class;
  return $self;
}

#Based on https://github.com/broadinstitute/gatk/blob/master/scripts/cnv_wdl/germline/cnv_germline_cohort_workflow.wdl
sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = get_parameter( $config, $section );

  my $java_option = $self->get_java_option( $config, $section, $memory );

  #parameter files
  my $gatk_singularity = get_param_file( $config->{$section}{gatk_singularity}, "gatk_singularity", 1 );

  my $intervals               = parse_param_file( $config, $section, "filtered_intervals",      1 );
  my $contig_ploidy_calls_dir = parse_param_file( $config, $section, "contig_ploidy_calls_dir", 1 );

  my $parameters = get_parameter_options(
    $config, $section, "--",
    [
      "p-alt",                                  "p-active",                         "cnv-coherence-length", "class-coherence-length", "max-copy-number", "max-bias-factors",
      "mapping-error-rate",                     "interval-psi-scale",               "sample-psi-scale",
      "depth-correction-tau",                   "log-mean-bias-standard-deviation", "init-ard-rel-unexplained-variance",
      "num-gc-bins",                            "gc-curve-standard-deviation",      "copy-number-posterior-expectation-mode",
      "enable-bias-factors",                    "active-class-padding-hybrid-mode", "learning-rate",
      "adamax-beta-1",                          "adamax-beta-2",                    "log-emission-samples-per-round",
      "log-emission-sampling-median-rel-error", "log-emission-sampling-rounds",
      "max-advi-iter-first-epoch",        "max-advi-iter-subsequent-epochs",     "min-training-epochs",
      "max-training-epochs",              "initial-temperature",                 "num-thermal-advi-iters",
      "convergence-snr-averaging-window", "convergence-snr-trigger-threshold",   "convergence-snr-countdown-window",
      "max-calling-iters",                "caller-update-convergence-threshold", "caller-internal-admixing-rate",
      "caller-external-admixing-rate",    "disable-annealing"
    ]
  );

  my $cohort_entity_id = $task_name;

  my $final_file = "${task_name}-contig-ploidy-calls.tar.gz";
  my $raw_files = get_raw_files( $config, $section );

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $inputOption = get_rawfiles_option( $raw_files, "--input" );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh "  
export HOME=$result_dir
export PYTHONPATH=

source activate gatk

cd $result_dir

gatk --java-options \"$java_option\" GermlineCNVCaller \\
  --run-mode COHORT \\
  -L ${intervals} \\
  $inputOption \\
  --contig-ploidy-calls $contig_ploidy_calls_dir \\
  --interval-merging-rule OVERLAPPING_ONLY \\
  --output . \\
  --output-prefix ${cohort_entity_id} \\
  --verbosity DEBUG $parameters \\
            
rm -rf .cache .conda .config .theano
";
  close($sh);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command );
  print $pbs "singularity run $gatk_singularity $shfile \n";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $result = { $task_name => filter_array( [ "${result_dir}/${task_name}-model", "${result_dir}/${task_name}-calls", "${result_dir}/${task_name}-tracking" ], $pattern ) };

  return $result;
}

1;
