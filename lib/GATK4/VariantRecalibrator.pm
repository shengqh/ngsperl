#!/usr/bin/perl
package GATK4::VariantRecalibrator;

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = get_parameter( $config, $section );

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

  my $vqsrMode = get_option( $config, $section, "vqsr_mode" );
  my $gvcf = get_option( $config, $section, "gvcf", 1 );

  if ($vqsrMode) {
    $gvcf   = 1;
    $hapmap_resource_vcf = get_param_file( $config->{$section}{hapmap_vcf}, "hapmap_vcf", 1 );
    $omni_resource_vcf   = get_param_file( $config->{$section}{omni_vcf}, "omni_vcf", 0 );
    $one_thousand_genomes_resource_vcf  = get_param_file( $config->{$section}{g1000_vcf}, "g1000_vcf", 0 );
    $mills_resource_vcf  = get_param_file( $config->{$section}{mills_vcf}, "mills_vcf", 1 );
    $dbsnp_resource_vcf  = get_param_file( $config->{$section}{dbsnp_vcf}, "dbsnp_vcf", 1 );
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

  my %vcfFiles = %{ get_raw_files( $config, $section ) };

  my $script = dirname(__FILE__) . "/fixLeftTrimDeletion.py";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $sites_only_variant_filtered_vcf       = $task_name . ".variant_filtered.sites_only.vcf.gz";
  
  my $indels_recalibration  = $task_name . ".indels.recal.vcf.gz";
  my $indels_tranches = $task_name . ".indels.tranches";
  
  my $snps_recalibration  = $task_name . ".snp.recal.vcf.gz";
  my $snps_tranches = $task_name . ".snp.tranches";
  
  my $indel_tranche_values = ["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"];
  my $snp_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0" ];
  my $recalibration_annotation_values = ["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"];

  my $indel_tranche_option = "-tranche " . join(" \\\n    -tranche ", @$indel_tranche_values);
  my $snp_tranche_option = "-tranche " . join(" \\\n    -tranche ", @$snp_tranche_values);
  my $recalibration_annotation_option = "-an " . join(" \\\n    -an ", @$recalibration_annotation_values);

  my $log_desc = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $indels_tranches, $init_command );
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
    @sample_names = sort keys %vcfFiles;
  }

  for my $sample_name ( @sample_names ) {
    my @sample_files = @{ $vcfFiles{$sample_name} };
    my $gvcfFile     = $sample_files[0];
    print $pbs "    -I $gvcfFile \\\n";
  }

  print $pbs "    -O $sites_only_variant_filtered_vcf

  tabix -p vcf $sites_only_variant_filtered_vcf
fi

if [[ -s $sites_only_variant_filtered_vcf && ! -s $indels_recalibration ]]; then
  echo IndelVariantRecalibrator=`date`
  gatk --java-options \"$java_option\" \\
    VariantRecalibrator \\
    -V ${sites_only_variant_filtered_vcf} \\
    -O $indels_recalibration \\
    --tranches-file $indels_tranches \\
    --trust-all-polymorphic \\
    $indel_tranche_option \\
    $recalibration_annotation_option \\
    --max-gaussians 4 \\
    $mills_option $axiomPoly_option -mode INDEL 
fi

if [[ -s $sites_only_variant_filtered_vcf && ! -s $snps_recalibration ]]; then
  echo SnpVariantRecalibrator=`date`
  gatk --java-options \"$java_option\" \\
    VariantRecalibrator \\
    -V ${sites_only_variant_filtered_vcf} \\
    -O $snps_recalibration \\
    --tranches-file $snps_tranches \\
    --trust-all-polymorphic \\
    $snp_tranche_option \\
    $recalibration_annotation_option \\
    --max-gaussians 6 \\
    ${hapmap_option} ${omni_option} ${g1000_option} ${dbsnp_option} -mode SNP
fi

rm -rf .conda

";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $indels_recalibration  = $result_dir . "/" . $task_name . ".indels.recal.vcf.gz";
  my $indels_tranches = $result_dir . "/" . $task_name . ".indels.tranches";
  
  my $snps_recalibration  = $result_dir . "/" . $task_name . ".snp.recal.vcf.gz";
  my $snps_tranches = $result_dir . "/" . $task_name . ".snp.tranches";

  my $result = {};
  $result->{$task_name} = filter_array( [$indels_recalibration, $indels_tranches, $snps_recalibration, $snps_tranches], $pattern );

  return $result;
}

1;
