#!/usr/bin/perl
package GATK4::VariantFilter;

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
  $self->{_suffix} = "_vf";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $vqsrMode = get_option( $config, $section, "vqsr_mode" );
  my $gvcf = get_option( $config, $section, "gvcf", 1 );
  
  my $indel_filter_level = get_option($config, $section, "indel_filter_level", 99.7);
  my $snp_filter_level = get_option($config, $section, "snp_filter_level", 99.7);

  #https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/joint-discovery-gatk4-local.wdl
  my $excess_het_threshold = 54.69;

  my $dbsnp_resource_vcf;
  my $mills_resource_vcf;
  my $axiomPoly_resource_vcf;
  my $hapmap_resource_vcf;
  my $omni_resource_vcf;
  my $one_thousand_genomes_resource_vcf;

  if ($vqsrMode) {
    my $species = get_option($config, $section, "species");
    my $isHuman = $species eq "homo_sapiens";
    $gvcf   = 1;
    $hapmap_resource_vcf = get_param_file( $config->{$section}{hapmap_vcf}, "hapmap_vcf", $isHuman );
    $omni_resource_vcf   = get_param_file( $config->{$section}{omni_vcf}, "omni_vcf", 0 );
    $one_thousand_genomes_resource_vcf  = get_param_file( $config->{$section}{g1000_vcf}, "g1000_vcf", 0 );
    $mills_resource_vcf  = get_param_file( $config->{$section}{mills_vcf}, "mills_vcf", $isHuman );
    $dbsnp_resource_vcf  = get_param_file( $config->{$section}{dbsnp_vcf}, "dbsnp_vcf", $isHuman );
    $axiomPoly_resource_vcf  = get_param_file( $config->{$section}{axiomPoly_vcf}, "axiomPoly_vcf", 0 );
  }
  
  my $hapmap_option = defined $hapmap_resource_vcf?"-resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap_resource_vcf} \\\n    ":"";
  my $omni_option = defined $omni_resource_vcf?"-resource:omni,known=false,training=true,truth=true,prior=12 ${omni_resource_vcf} \\\n    ":"";
  my $g1000_option = defined $one_thousand_genomes_resource_vcf? "-resource:1000G,known=false,training=true,truth=false,prior=10 ${one_thousand_genomes_resource_vcf} \\\n    ":"";
  my $axiomPoly_option = defined $axiomPoly_resource_vcf? "-resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiomPoly_resource_vcf} \\\n    ":"";
  my $mills_option = defined $mills_resource_vcf? "-resource:mills,known=false,training=true,truth=true,prior=12 ${mills_resource_vcf} \\\n    ":"";
  my $dbsnp_option = defined $dbsnp_resource_vcf? "-resource:dbsnp,known=true,training=false,truth=false,prior=2 ${dbsnp_resource_vcf} \\\n    ":"";

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my $java_option = $self->get_java_option($config, $section, $memory);
  $self->get_docker_value(1);

  my $vcfFiles = get_raw_files( $config, $section );

  my $script = dirname(__FILE__) . "/fixLeftTrimDeletion.py";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $indel_tranche_values = ["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"];
  my $snp_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0" ];
  my $recalibration_annotation_values = ["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"];

  my $indel_tranche_option = "-tranche " . join(" \\\n    -tranche ", @$indel_tranche_values);
  my $snp_tranche_option = "-tranche " . join(" \\\n    -tranche ", @$snp_tranche_values);
  my $recalibration_annotation_option = "-an " . join(" \\\n    -an ", @$recalibration_annotation_values);
  
  my $indel_recalibration_tmp_vcf = $task_name . ".indels.recal.tmp.vcf.gz";
  my $recalibrated_vcf_filename = $task_name . ".indels.snp.recal.vcf.gz";

  my $pass_file = $task_name . ".indels.snp.recal.pass.vcf.gz";
  my $split_file = $task_name . ".indels.snp.recal.pass.split.vcf";
  my $left_trim_file = $task_name . ".indels.snp.recal.pass.norm.vcf";
  my $fix_file = $task_name . ".indels.snp.recal.pass.norm.nospan.vcf";
  my $final_file = $task_name . ".indels.snp.recal.pass.norm.nospan.vcf.gz";

  my ($pbs_file, $pbs, $variant_filtered_vcf, $sites_only_variant_filtered_vcf, $rmlist) = start_variant_filter_pbs($self, $task_name, $pbs_dir, $pbs_desc, $log_dir, $thread, $cluster, $path_file, $result_dir, $init_command, $java_option, $faFile, $vcfFiles, $dbsnp_resource_vcf, $excess_het_threshold, undef, $final_file);

  my ($indels_recalibration, $indels_tranches, $snps_recalibration, $snps_tranches, $recal_rmlist) = add_indel_snv_recalibrator_pbs($self, $config, $section, $pbs, $task_name, $sites_only_variant_filtered_vcf, $memory);

  my ($apply_vqsr_rmlist) = add_apply_vqsr_pbs($self, $config, $section, $pbs, $task_name, $variant_filtered_vcf, $indels_recalibration, $indels_tranches, $snps_recalibration, $snps_tranches, $pass_file, $option, $thread, $memory );

  my ($left_trim_rmlist) = add_left_trim_pbs($self, $config, $section, $pbs, $task_name, $pass_file, $final_file);

  print $pbs "
if [[ -s $final_file ]]; then
  rm $rmlist \\
    $recal_rmlist \\
    $apply_vqsr_rmlist \\
    $indels_recalibration ${indels_recalibration}.tbi $indels_tranches  \\
    $snps_recalibration ${snps_recalibration}.tbi $snps_tranches \\
    $left_trim_rmlist \\
    .conda
fi

";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $final_file = $task_name . ".indels.snp.recal.pass.norm.nospan.vcf.gz";
  my @result_files = ();
  push( @result_files, $result_dir . "/" . $final_file );

  my $result = {};
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
