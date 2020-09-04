#!/usr/bin/perl
package GATK4::VariantFilterRecalibrator;

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
use GATK4::GATK4UniqueTask;
use GATK4::VariantFilterUtils;

our @ISA = qw(GATK4::GATK4UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_vr";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $vcfFiles = get_raw_files($config, $section);
  my $java_option = $self->get_java_option($config, $section, $memory);
  $self->get_docker_value(1);

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $sites_only_variant_filtered_vcf       = $task_name . ".variant_filtered.sites_only.vcf.gz";
  my $snps_tranches = $task_name . ".snp.tranches";
  
  my $log_desc = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $snps_tranches, $init_command );
  print $pbs "  
if [ ! -s $sites_only_variant_filtered_vcf ]; then
  echo GatherVcfsCloud=`date` 
  gatk --java-options \"$java_option\" \\
    GatherVcfsCloud \\
    --ignore-safety-checks \\
    --gather-type BLOCK \\
";

  my @sample_names = ();
  if (has_option($config, $section, "chromosome_names")){
    my $chromosomeStr = get_option($config, $section, "chromosome_names");
    my @chromosomes = split /,/, $chromosomeStr;
    for my $chr (@chromosomes) {
      my $chrTaskName = $task_name . "." . $chr;
      push @sample_names, $chrTaskName;
    }
  }else{
    @sample_names = sort keys %$vcfFiles;
  }

  for my $sample_name ( @sample_names ) {
    my @sample_files = @{ $vcfFiles->{$sample_name} };
    my $gvcfFile     = $sample_files[0];
    print $pbs "    -I $gvcfFile \\\n";
  }

  print $pbs "    -O $sites_only_variant_filtered_vcf

  tabix -p vcf $sites_only_variant_filtered_vcf
fi
";

  add_indel_snv_recalibrator_pbs($self, $config, $section, $pbs, $task_name, $sites_only_variant_filtered_vcf, $memory);

  print $pbs "
rm -rf .conda

";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $indels_recalibration  = $result_dir . "/" . $task_name . ".indels.recal.vcf.gz";
  my $indels_tranches = $result_dir . "/" . $task_name . ".indels.tranches";
  
  my $snps_recalibration  = $result_dir . "/" . $task_name . ".snp.recal.vcf.gz";
  my $snps_tranches = $result_dir . "/" . $task_name . ".snp.tranches";

  my $result = {};
  $result->{$task_name} = filter_array( [$indels_recalibration, $indels_tranches, $snps_recalibration, $snps_tranches], $pattern );

  return $result;
}

1;
