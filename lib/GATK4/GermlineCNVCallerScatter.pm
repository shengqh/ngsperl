#!/usr/bin/perl
package GATK4::GermlineCNVCallerScatter;

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
use GATK4::GATK4Task;

our @ISA = qw(GATK4::GATK4Task);

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $java_option = $self->get_java_option( $config, $section, $memory );

  #parameter files
  $self->get_docker_value(1);

  my $contig_ploidy_calls_dir = parse_param_file( $config, $section, "contig_ploidy_calls_dir", 1 );

  my $parameters = get_parameter_options(
    $config, $section, "--",
    [
      "p-alt",                            "p-active",                               "cnv-coherence-length",                   "class-coherence-length",               #
      "max-copy-number",                  "max-bias-factors",                       "mapping-error-rate",                     "interval-psi-scale",                   #
      "sample-psi-scale",                 "depth-correction-tau",                   "log-mean-bias-standard-deviation",       "init-ard-rel-unexplained-variance",    #
      "num-gc-bins",                      "gc-curve-standard-deviation",            "copy-number-posterior-expectation-mode", "enable-bias-factors",                  #
      "active-class-padding-hybrid-mode", "learning-rate",                          "adamax-beta-1",                          "adamax-beta-2",                        #
      "log-emission-samples-per-round",   "log-emission-sampling-median-rel-error", "log-emission-sampling-rounds",           "max-advi-iter-first-epoch",            #
      "max-advi-iter-subsequent-epochs",  "min-training-epochs",                    "max-training-epochs",                    "initial-temperature",                  #
      "num-thermal-advi-iters",           "convergence-snr-averaging-window",       "convergence-snr-trigger-threshold",      "convergence-snr-countdown-window",     #
      "max-calling-iters", "caller-update-convergence-threshold", "caller-internal-admixing-rate",                                                                    #
      "caller-external-admixing-rate", "disable-annealing"
    ],
    [
      "1e-6",    "1e-2",  "10000.0",                                                                                                                                  #
      "10000.0", "5",     "5",                                                                                                                                        #
      "0.01",    "0.001", "0.0001",                                                                                                                                   #
      "10000.0", "0.1",   "0.1",                                                                                                                                      #
      "20",      "1.0",   "HYBRID",                                                                                                                                   #
      "true",    "50000", "0.05",                                                                                                                                     #
      "0.9",     "0.99",  "50",                                                                                                                                       #
      "0.005", "10",                                                                                                                                                  #
      "5000",  "100",                                                                                                                                                 #
      "10",    "100",                                                                                                                                                 #
      "2.0",   "2500",                                                                                                                                                #
      "500",   "0.1", "10",                                                                                                                                           #
      "10",    "0.001", "0.75",                                                                                                                                       #
      "1.00",  "false"
    ]
  );

  my $intervals = get_raw_files( $config, $section);
  my $count_files =  get_raw_files( $config, $section, "counts");
  my $inputOption = get_rawfiles_option( $count_files, "--input", " \\\n  " );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $scatterName (sort keys %$intervals){
    my $cur_dir  = create_directory_or_die( $result_dir . "/$scatterName" );
    my $final_file = "$cur_dir/gcc-calls/SAMPLE_0/sample_name.txt";

    my $cohort_entity_id = $task_name . "_" . $scatterName;
  
    my $interval_file = $intervals->{$scatterName}[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $scatterName );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $scatterName );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "if [[ ! -s $final_file ]]; then
  \$MYCMD ./$pbs_name 
fi
";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file, $init_command );
    print $pbs "  

cd $cur_dir

export MKL_NUM_THREADS=$thread
export OMP_NUM_THREADS=$thread

gatk --java-options \"$java_option\" GermlineCNVCaller $option \\
  --run-mode COHORT \\
  -L ${interval_file} $inputOption \\
  --contig-ploidy-calls $contig_ploidy_calls_dir \\
  --interval-merging-rule OVERLAPPING_ONLY \\
  --output . \\
  --output-prefix gcc \\
  --verbosity DEBUG $parameters
            
rm -rf .cache .conda .config .theano
";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks. \n ";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $result = {};
  my $intervals = get_raw_files( $config, $section);
  for my $scatterName (sort keys %$intervals) {
    $result->{$scatterName} = filter_array( [ "${result_dir}/${scatterName}/gcc-model", "${result_dir}/${scatterName}/gcc-calls", "${result_dir}/${scatterName}/gcc-tracking" ], $pattern );
  }
  return $result;
}

1;
