#!/usr/bin/perl
package GATK4::VariantFilterScatterChromosome;

use strict;
use warnings;
use File::Basename;
use List::Util qw[min];
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;
use GATK4::GATK4ChromosomeTask;
use GATK4::VariantFilterUtils;

our @ISA = qw(GATK4::GATK4ChromosomeTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_vfc";
  $self->{_depend_all} = 1;
  bless $self, $class;
  return $self;
}

#get input file list
sub get_sample_names {
  my ($self, $config, $section) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  return [$task_name];  
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $vqsrMode = get_option( $config, $section, "vqsr_mode" );
  my $gvcf = get_option( $config, $section, "gvcf", 1 );

  my $chromosomeStr = get_option($config, $section, "chromosome_names");
  my @chromosomes = split /,/, $chromosomeStr;
  
  my $indel_filter_level = get_option($config, $section, "indel_filter_level", 99.7);
  my $snp_filter_level = get_option($config, $section, "snp_filter_level", 99.7);

  #https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/joint-discovery-gatk4-local.wdl
  my $excess_het_threshold = 54.69;

  my $dbsnp_resource_vcf;

  if ($vqsrMode) {
    $gvcf   = 1;
    $dbsnp_resource_vcf  = get_param_file( $config->{$section}{dbsnp_vcf}, "dbsnp_vcf", 1 );
  }
  
  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my $java_option = $self->get_java_option($config, $section, $memory);
  $self->get_docker_value(1);

  my $vcfFiles = get_raw_files( $config, $section );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $chr (@chromosomes) {
    my ($pbs_file, $pbs, $variant_filtered_vcf, $final_file, $rmlist) = start_variant_filter_pbs($self, $task_name, $pbs_dir, $pbs_desc, $log_dir, $thread, $cluster, $path_file, $result_dir, $init_command, $java_option, $faFile, $vcfFiles, $dbsnp_resource_vcf, $excess_het_threshold, $chr);
    my $pbs_name = basename($pbs_file);
    print $sh "\$MYCMD ./$pbs_name \n";

    print $pbs "
if [[ -s $final_file ]]; then
  rm $rmlist .conda
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
  my $variant_filtered_vcf = $result_dir . "/" . $key_name . ".variant_filtered.vcf.gz";
  my $sites_only_variant_filtered_vcf = $result_dir . "/" . $key_name . ".variant_filtered.sites_only.vcf.gz";
  return([$variant_filtered_vcf, $sites_only_variant_filtered_vcf]);
}

1;
