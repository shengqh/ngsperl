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
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = get_parameter( $config, $section );

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

  my $mergedFile                    = $task_name . ".merged.vcf.gz";
  my $callFile                      = $task_name . ".raw.vcf.gz";
  my $variant_filtered_vcf = $task_name . ".variant_filtered.vcf.gz";
  my $sites_only_variant_filtered_vcf       = $task_name . ".sites_only.variant_filtered.vcf.gz";
  
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
  
  my $indel_recalibration_tmp_vcf = $task_name . ".indels.recal.tmp.vcf.gz";
  my $recalibrated_vcf_filename = $task_name . ".indels.snp.recal.vcf.gz";

  my $pass_file = $task_name . ".indels.snp.recal.pass.vcf.gz";
  my $split_file = $task_name . ".indels.snp.recal.pass.split.vcf";
  my $left_trim_file = $task_name . ".indels.snp.recal.pass.norm.vcf";
  my $fix_file = $task_name . ".indels.snp.recal.pass.norm.nospan.vcf";
  my $final_file = $task_name . ".indels.snp.recal.pass.norm.nospan.vcf.gz";

  my $reader_threads = min( 5, $thread );

  my $log_desc = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command );
  print $pbs "  
cd $result_dir

if [ ! -s $mergedFile ]; then
  echo CombineGVCFs=`date` 
  gatk --java-options \"$java_option\" \\
    CombineGVCFs \\
    -R $faFile \\
";

  for my $sample_name ( sort keys %vcfFiles ) {
    my @sample_files = @{ $vcfFiles{$sample_name} };
    my $gvcfFile     = $sample_files[0];
    print $pbs "    -V $gvcfFile \\\n";
  }

  print $pbs "    -O $mergedFile
fi

if [[ -s $mergedFile && ! -s $callFile ]]; then
  echo GenotypeGVCFs=`date` 
  gatk --java-options \"$java_option\" \\
    GenotypeGVCFs \\
    -R $faFile \\
    -D $dbsnp_resource_vcf \\
    -O $callFile \\
    -V $mergedFile 
fi

if [[ -s $callFile && ! -s $variant_filtered_vcf ]]; then
  echo VariantFiltration=`date` 
  gatk --java-options \"$java_option\" \\
    VariantFiltration \\
    --filter-expression \"ExcessHet > ${excess_het_threshold}\" \\
    --filter-name ExcessHet \\
    -O $variant_filtered_vcf \\
    -V ${callFile}
fi

if [[ -s $variant_filtered_vcf && ! -s $sites_only_variant_filtered_vcf ]]; then
  echo MakeSitesOnlyVcf=`date` 
  gatk --java-options \"$java_option\" \\
    MakeSitesOnlyVcf \\
    --INPUT $variant_filtered_vcf \\
    --OUTPUT $sites_only_variant_filtered_vcf
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

if [[ -s $indels_recalibration && ! -s $indel_recalibration_tmp_vcf ]]; then
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

if [[ -s $pass_file && ! -s $left_trim_file ]]; then
  echo LeftAlignAndNorm=`date`
  bcftools norm -m- -o $split_file $pass_file 
  bcftools norm -f $faFile -o $left_trim_file $split_file 
fi

if [[ -s $left_trim_file && ! -s $final_file ]]; then
  echo noSpanDeletion=`date`
  python $script -i $left_trim_file -o $fix_file
  bgzip $fix_file
  tabix -p vcf $final_file
fi

if [[ -s $final_file ]]; then
  rm $mergedFile ${mergedFile}.tbi \\
    $callFile ${callFile}.tbi \\
    $variant_filtered_vcf ${variant_filtered_vcf}.tbi \\
    $sites_only_variant_filtered_vcf ${sites_only_variant_filtered_vcf}.tbi \\
    $indels_recalibration ${indels_recalibration}.tbi \\
    $snps_recalibration ${snps_recalibration}.tbi \\
    $indel_recalibration_tmp_vcf ${indel_recalibration_tmp_vcf}.tbi \\
    $split_file ${split_file}.idx $left_trim_file ${left_trim_file}.idx
fi

";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $final_file = $task_name . ".indels.snp.recal.pass.norm.nospan.vcf.gz";
  my @result_files = ();
  push( @result_files, $result_dir . "/" . $final_file );

  my $result = {};
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
