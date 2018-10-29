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

  my $java_option = $self->get_java_option($config, $section, $memory);

  #parameter files
  my $gatk_singularity = get_param_file( $config->{$section}{gatk_singularity}, "gatk_singularity", 1 );

  my $intervals            = parse_param_file( $config, $section, "filtered_intervals", 1);
  my $contig_ploidy_calls_dir            = parse_param_file( $config, $section, "contig_ploidy_calls_dir", 1);

  my $p_alt                  = get_option( $config, $section, "p_alt",                  1e-6 );
  my $p_active               = get_option( $config, $section, "p_active",               1e-2 );
  my $cnv_coherence_length   = get_option( $config, $section, "cnv_coherence_length",   10000.0 );
  my $class_coherence_length = get_option( $config, $section, "class_coherence_length", 10000.0 );

  my $max_copy_number                        = get_option( $config, $section, "max_copy_number",                        5 );
  my $max_bias_factors                       = get_option( $config, $section, "max_bias_factors",                       5 );
  my $mapping_error_rate                     = get_option( $config, $section, "mapping_error_rate",                     0.01 );
  my $interval_psi_scale                     = get_option( $config, $section, "interval_psi_scale",                     0.001 );
  my $sample_psi_scale                       = get_option( $config, $section, "sample_psi_scale",                       0.0001 );
  my $depth_correction_tau                   = get_option( $config, $section, "depth_correction_tau",                   10000.0 );
  my $log_mean_bias_standard_deviation       = get_option( $config, $section, "log_mean_bias_standard_deviation",       0.1 );
  my $init_ard_rel_unexplained_variance      = get_option( $config, $section, "init_ard_rel_unexplained_variance",      0.1 );
  my $num_gc_bins                            = get_option( $config, $section, "num_gc_bins",                            20 );
  my $gc_curve_standard_deviation            = get_option( $config, $section, "gc_curve_standard_deviation",            1.0 );
  my $copy_number_posterior_expectation_mode = get_option( $config, $section, "copy_number_posterior_expectation_mode", "HYBRID" );

  my $enable_bias_factors                    = get_option( $config, $section, "enable_bias_factors",                    "true" );
  my $active_class_padding_hybrid_mode       = get_option( $config, $section, "active_class_padding_hybrid_mode",       50000 );
  my $learning_rate                          = get_option( $config, $section, "learning_rate",                          0.05 );
  my $adamax_beta_1                          = get_option( $config, $section, "adamax_beta_1",                          0.9 );
  my $adamax_beta_2                          = get_option( $config, $section, "adamax_beta_2",                          0.99 );
  my $log_emission_samples_per_round         = get_option( $config, $section, "log_emission_samples_per_round",         50 );
  my $log_emission_sampling_median_rel_error = get_option( $config, $section, "log_emission_sampling_median_rel_error", 0.005 );
  my $log_emission_sampling_rounds           = get_option( $config, $section, "log_emission_sampling_rounds",           10 );
  my $max_advi_iter_first_epoch              = get_option( $config, $section, "max_advi_iter_first_epoch",              5000 );
  my $max_advi_iter_subsequent_epochs        = get_option( $config, $section, "max_advi_iter_subsequent_epochs",        100 );
  my $min_training_epochs                    = get_option( $config, $section, "min_training_epochs",                    10 );
  my $max_training_epochs                    = get_option( $config, $section, "max_training_epochs",                    100 );
  my $initial_temperature                    = get_option( $config, $section, "initial_temperature",                    2.0 );
  my $num_thermal_advi_iters                 = get_option( $config, $section, "num_thermal_advi_iters",                 2500 );
  my $convergence_snr_averaging_window       = get_option( $config, $section, "convergence_snr_averaging_window",       500 );
  my $convergence_snr_trigger_threshold      = get_option( $config, $section, "convergence_snr_trigger_threshold",      0.1 );
  my $convergence_snr_countdown_window       = get_option( $config, $section, "convergence_snr_countdown_window",       10 );
  my $max_calling_iters                      = get_option( $config, $section, "max_calling_iters",                      10 );
  my $caller_update_convergence_threshold    = get_option( $config, $section, "caller_update_convergence_threshold",    0.001 );
  my $caller_internal_admixing_rate          = get_option( $config, $section, "caller_internal_admixing_rate",          0.75 );
  my $caller_external_admixing_rate          = get_option( $config, $section, "caller_external_admixing_rate",          1.00 );
  my $disable_annealing                      = get_option( $config, $section, "disable_annealing",                      "false" );

  #my $ = get_option( $config, $section, "",  );

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
            --verbosity DEBUG \\
            --p-alt ${p_alt} \\
            --p-active ${p_active} \\
            --cnv-coherence-length ${cnv_coherence_length} \\
            --class-coherence-length ${class_coherence_length} \\
            --max-copy-number ${max_copy_number} \\
            --max-bias-factors ${max_bias_factors} \\
            --mapping-error-rate ${mapping_error_rate} \\
            --interval-psi-scale ${interval_psi_scale} \\
            --sample-psi-scale ${sample_psi_scale} \\
            --depth-correction-tau ${depth_correction_tau} \\
            --log-mean-bias-standard-deviation ${log_mean_bias_standard_deviation} \\
            --init-ard-rel-unexplained-variance ${init_ard_rel_unexplained_variance} \\
            --num-gc-bins ${num_gc_bins} \\
            --gc-curve-standard-deviation ${gc_curve_standard_deviation} \\
            --copy-number-posterior-expectation-mode ${copy_number_posterior_expectation_mode} \\
            --enable-bias-factors ${enable_bias_factors} \\
            --active-class-padding-hybrid-mode ${active_class_padding_hybrid_mode} \\
            --learning-rate ${learning_rate} \\
            --adamax-beta-1 ${adamax_beta_1} \\
            --adamax-beta-2 ${adamax_beta_2} \\
            --log-emission-samples-per-round ${log_emission_samples_per_round} \\
            --log-emission-sampling-median-rel-error ${log_emission_sampling_median_rel_error} \\
            --log-emission-sampling-rounds ${log_emission_sampling_rounds} \\
            --max-advi-iter-first-epoch ${max_advi_iter_first_epoch} \\
            --max-advi-iter-subsequent-epochs ${max_advi_iter_subsequent_epochs} \\
            --min-training-epochs ${min_training_epochs} \\
            --max-training-epochs ${max_training_epochs} \\
            --initial-temperature ${initial_temperature} \\
            --num-thermal-advi-iters ${num_thermal_advi_iters} \\
            --convergence-snr-averaging-window ${convergence_snr_averaging_window} \\
            --convergence-snr-trigger-threshold ${convergence_snr_trigger_threshold} \\
            --convergence-snr-countdown-window ${convergence_snr_countdown_window} \\
            --max-calling-iters ${max_calling_iters} \\
            --caller-update-convergence-threshold ${caller_update_convergence_threshold} \\
            --caller-internal-admixing-rate ${caller_internal_admixing_rate} \\
            --caller-external-admixing-rate ${caller_external_admixing_rate} \\
            --disable-annealing ${disable_annealing}
            
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

  my $result = {
    $task_name =>
      filter_array( [ "${result_dir}/${task_name}-model", "${result_dir}/${task_name}-calls", "${result_dir}/${task_name}-tracking" ], $pattern )
  };

  return $result;
}

1;
