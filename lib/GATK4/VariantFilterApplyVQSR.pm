#!/usr/bin/perl
package GATK4::VariantFilterApplyVQSR;

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
use CQS::TaskUtils;
use GATK4::GATK4ChromosomeTask;
use GATK4::VariantFilterUtils;

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
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  return [$task_name];  
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $java_option = $self->get_java_option($config, $section, $memory);
  $self->get_docker_value(1);

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

  for my $chrTaskName (sort keys %$vcfFiles) {
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $chrTaskName );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $chrTaskName );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";

    my $variant_filtered_vcf = $vcfFiles->{$chrTaskName}[0];
    my $pass_file = $chrTaskName . ".indels.snp.recal.pass.vcf.gz";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $pass_file, $init_command );

    my ($rmlist) = add_apply_vqsr_pbs($self, $config, $section, $pbs, $chrTaskName, $variant_filtered_vcf, 
      $indels_recalibration, $indels_tranches, $snps_recalibration, $snps_tranches, $pass_file, $option, $thread, $memory );
print $pbs "
if [[ -s $pass_file ]]; then
  rm -rf $rmlist \\
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
