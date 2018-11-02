#!/usr/bin/perl
package GATK4::PostprocessGermlineCNVCalls;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_pgcc";
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

  my $contig_ploidy_calls_dir = parse_param_file( $config, $section, "contig_ploidy_calls_dir", 1 );

  my $calls_shard_path = get_raw_files( $config, $section, "calls_shard_path" );
  my $calls_args = get_rawfiles_option( $calls_shard_path, "--calls-shard-path" );

  my $model_shard_path = get_raw_files( $config, $section, "model_shard_path" );
  my $model_args = get_rawfiles_option( $model_shard_path, "--model-shard-path" );
  
  my $parameters = get_parameter_options($config, $section, "--", "autosomal-ref-copy-number");

  #make PBS
  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  my @sample_names = sort keys %raw_files;
  for my $i ( 0 .. $#sample_names ) {
    my $sample_name = $sample_names[$i];

    my $sample_dir = create_directory_or_die($result_dir . "/" . $sample_name );
    
    my $genotyped_intervals_vcf_filename = $sample_name . ".genotyped_intervals.vcf.gz";
    my $genotyped_segments_vcf_filename  = $sample_name . ".genotyped_segments.vcf.gz";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $shsamplefile = $self->get_task_filename( $pbs_dir, $sample_name );
    open( my $shsample, ">$shsamplefile" ) or die "Cannot create $shsamplefile";
    print $shsample "  
export HOME=$result_dir
export PYTHONPATH=

source activate gatk

cd $sample_dir

gatk --java-options \"$java_option\" PostprocessGermlineCNVCalls \\
  $calls_args \\
  $model_args \\
  --sample-index $i \\
  --allosomal-contig X \\
  --allosomal-contig Y $parameters \\
  --contig-ploidy-calls $contig_ploidy_calls_dir \\
  --output-genotyped-intervals $genotyped_intervals_vcf_filename \\
  --output-genotyped-segments $genotyped_segments_vcf_filename
            
rm -rf .cache .conda .config .theano

";
    close($shsample);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $genotyped_intervals_vcf_filename, $init_command );
    print $pbs "singularity run $gatk_singularity $shsamplefile \n";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all GATK CNV tasks. \n ";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $sample_dir = create_directory_or_die($result_dir . "/" . $sample_name );
    my $genotyped_intervals_vcf_filename = $sample_name . ".genotyped_intervals.vcf.gz";
    my $genotyped_segments_vcf_filename  = $sample_name . ".genotyped_segments.vcf.gz";
    $result->{$sample_name} = filter_array( [ "${sample_dir}/$genotyped_intervals_vcf_filename", "${sample_dir}/$genotyped_segments_vcf_filename" ], $pattern );
  }

  return $result;
}

1;
