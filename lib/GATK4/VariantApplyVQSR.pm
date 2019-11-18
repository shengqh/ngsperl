#!/usr/bin/perl
package GATK4::VariantApplyVQSR;

use strict;
use warnings;
use File::Basename;
use List::Util qw[min];
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use GATK4::GATK4ChromosomeTask;

our @ISA = qw(GATK4::GATK4ChromosomeTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_vfa";
  bless $self, $class;
  return $self;
}

sub get_sample_names {
  my ($self, $config, $section) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );
  return [$task_name];  
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = get_parameter( $config, $section );

  my $java_option = $self->get_java_option($config, $section, $memory);
  $self->get_docker_value(1);

  my $chromosomeStr = get_option($config, $section, "chromosome_names");
  my @chromosomes = split /,/, $chromosomeStr;

  my $vcfFiles = get_raw_files( $config, $section );

  my $indels_recalibration  = parse_param_file($config, $section, "indels_recalibration");
  my $indels_tranches = parse_param_file($config, $section, "indels_tranches");
  my $snps_recalibration  = parse_param_file($config, $section, "snps_recalibration");
  my $snps_tranches = parse_param_file($config, $section, "snps_tranches");
  
  my $indel_filter_level = get_option($config, $section, "indel_filter_level", 99.7);
  my $snp_filter_level = get_option($config, $section, "snp_filter_level", 99.7);

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $chr (@chromosomes) {
    my $chrTaskName = $self->get_key_name($task_name, $chr);

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $chrTaskName );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $chrTaskName );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $variant_filtered_vcf = $vcfFiles->{$chrTaskName}[0];
    my $reader_threads = min( 5, $thread );
    my $log_desc = $cluster->get_log_description($log);

    my $indel_recalibration_tmp_vcf = $chrTaskName . ".indels.recal.tmp.vcf.gz";
    my $recalibrated_vcf_filename = $chrTaskName . ".indels.snp.recal.vcf.gz";
    my $pass_file = $chrTaskName . ".indels.snp.recal.pass.vcf.gz";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $pass_file, $init_command );

    print $pbs "  

if [[ ! -s $indel_recalibration_tmp_vcf ]]; then
  echo IndelApplyVQSR=`date`
  gatk --java-options \"$java_option\" \\
    ApplyVQSR \\
    -O $indel_recalibration_tmp_vcf \\
    -V $variant_filtered_vcf \\
    --recal-file $indels_recalibration \\
    --tranches-file $indels_tranches \\
    --truth-sensitivity-filter-level $indel_filter_level \\
    --create-output-variant-index true \\
    -mode INDEL
fi

if [[ -s $snps_recalibration && ! -s $recalibrated_vcf_filename ]]; then
  echo SnpApplyVQSR=`date`
  gatk --java-options \"$java_option\" \\
    ApplyVQSR \\
    -O $recalibrated_vcf_filename \\
    -V $indel_recalibration_tmp_vcf \\
    --recal-file $snps_recalibration \\
    --tranches-file $snps_tranches \\
    --truth-sensitivity-filter-level $snp_filter_level \\
    --create-output-variant-index true \\
    -mode SNP
fi

if [[ -s $recalibrated_vcf_filename && ! -s $pass_file ]]; then
  echo SelectVariant=`date`
  gatk --java-options \"$java_option\" \\
    SelectVariants \\
    -O $pass_file \\
    -V $recalibrated_vcf_filename \\
    --exclude-filtered
fi

if [[ -s $pass_file ]]; then
  rm -rf $indel_recalibration_tmp_vcf ${indel_recalibration_tmp_vcf}.tbi \\
    .conda
fi

";
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks. \n ";
}

sub get_result_files {
  my ( $self, $config, $section, $result_dir, $sample_name, $scatter_name, $key_name ) = @_;
  my $pass_file = $result_dir . "/" . $key_name . ".indels.snp.recal.pass.vcf.gz";
  return([$pass_file]);
}

1;
