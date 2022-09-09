#!/usr/bin/perl
package GATK4::VariantFilterUtils;

use strict;
use warnings;
use File::Basename;
use List::Util qw[min];
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw( start_variant_filter_pbs add_indel_snv_recalibrator_pbs add_apply_vqsr_pbs add_left_trim_pbs)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

sub start_variant_filter_pbs {
  my ($classobj, $task_name, $pbs_dir, $pbs_desc, $log_dir, $thread, $cluster, $path_file, $result_dir, $init_command, $java_option, $faFile, $vcfFiles, $dbsnp_resource_vcf, $excess_het_threshold, $chr) = @_;

  my $chrTaskName = (defined $chr) ? get_key_name($task_name, $chr) : $task_name;
  my $chrOption = (defined $chr) ? " -L $chr " : "";

  my $pbs_file = $classobj->get_pbs_filename( $pbs_dir, $chrTaskName );
  my $log      = $classobj->get_log_filename( $log_dir, $chrTaskName );

  my $mergedFile                    = $chrTaskName . ".merged.vcf.gz";
  my $callFile                      = $chrTaskName . ".raw.vcf.gz";
  my $variant_filtered_vcf = $chrTaskName . ".variant_filtered.vcf.gz";
  my $final_file       = $chrTaskName . ".variant_filtered.sites_only.vcf.gz";
  
  my $reader_threads = min( 5, $thread );

  my $log_desc = $cluster->get_log_description($log);

  my $pbs = $classobj->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command );
  print $pbs "  
if [ ! -s $mergedFile ]; then
  echo CombineGVCFs=`date` 
  gatk --java-options \"$java_option\" \\
    CombineGVCFs \\
    -R $faFile $chrOption\\
";

  for my $sample_name ( sort keys %$vcfFiles ) {
    my @sample_files = @{ $vcfFiles->{$sample_name} };
    my $gvcfFile     = $sample_files[0];
    print $pbs "      -V $gvcfFile \\\n";
  }

  print $pbs "      -O $mergedFile

  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $chrTaskName.CombinedGVCFs.failed
    rm $mergedFile
  fi
fi

if [[ -s $mergedFile && ! -s $callFile ]]; then
  echo GenotypeGVCFs=`date` 
  gatk --java-options \"$java_option\" \\
    GenotypeGVCFs \\
    -R $faFile \\
    -D $dbsnp_resource_vcf \\
    -O $callFile \\
    -V $mergedFile 

  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $chrTaskName.GenotypeGVCFs.failed
    rm $callFile
  fi
fi

if [[ -s $callFile && ! -s $variant_filtered_vcf ]]; then
  echo VariantFiltration=`date` 
  gatk --java-options \"$java_option\" \\
    VariantFiltration \\
    --filter-expression \"ExcessHet > ${excess_het_threshold}\" \\
    --filter-name ExcessHet \\
    -O $variant_filtered_vcf \\
    -V ${callFile}

  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $chrTaskName.VariantFiltration.failed
    rm $variant_filtered_vcf
  fi
fi

if [[ -s $variant_filtered_vcf && ! -s $final_file ]]; then
  echo MakeSitesOnlyVcf=`date` 
  gatk --java-options \"$java_option\" \\
    MakeSitesOnlyVcf \\
    --INPUT $variant_filtered_vcf \\
    --OUTPUT $final_file

  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $chrTaskName.MakeSitesOnlyVcf.failed
    rm $final_file
  fi
fi
";

  my $rmlist = "$mergedFile ${mergedFile}.tbi \\\n    $callFile ${callFile}.tbi";
  return ($pbs_file, $pbs, $variant_filtered_vcf, $final_file, $rmlist);
}

sub add_indel_snv_recalibrator_pbs {
  my ( $classobj, $config, $section, $pbs, $task_name, $sites_only_variant_filtered_vcf, $memory ) = @_;

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

  #my $downsampleFactor = get_option($config, $section, "SNP_VQSR_downsampleFactor", 10);

  my $vqsrMode = get_option( $config, $section, "vqsr_mode" );
  my $gvcf = get_option( $config, $section, "gvcf", 1 );

  if ($vqsrMode) {
    $gvcf   = 1;
    $hapmap_resource_vcf = get_param_file( $config->{$section}{hapmap_vcf}, "hapmap_vcf", 0 );
    $omni_resource_vcf   = get_param_file( $config->{$section}{omni_vcf}, "omni_vcf", 0 );
    $one_thousand_genomes_resource_vcf  = get_param_file( $config->{$section}{g1000_vcf}, "g1000_vcf", 0 );
    $mills_resource_vcf  = get_param_file( $config->{$section}{mills_vcf}, "mills_vcf", 0 );
    $dbsnp_resource_vcf  = get_param_file( $config->{$section}{dbsnp_vcf}, "dbsnp_vcf", 1 );
    $axiomPoly_resource_vcf  = get_param_file( $config->{$section}{axiomPoly_vcf}, "axiomPoly_vcf", 0 );
  }
  
  my $hapmap_option = defined $hapmap_resource_vcf?"--resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap_resource_vcf} \\\n    ":"";
  my $omni_option = defined $omni_resource_vcf?"--resource:omni,known=false,training=true,truth=true,prior=12 ${omni_resource_vcf} \\\n    ":"";
  my $g1000_option = defined $one_thousand_genomes_resource_vcf? "--resource:1000G,known=false,training=true,truth=false,prior=10 ${one_thousand_genomes_resource_vcf} \\\n    ":"";
  my $axiomPoly_option = defined $axiomPoly_resource_vcf? "--resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiomPoly_resource_vcf} \\\n    ":"";
  my $mills_option = defined $mills_resource_vcf? "--resource:mills,known=false,training=true,truth=true,prior=12 ${mills_resource_vcf} \\\n    ":"";
  my $dbsnp_option = defined $dbsnp_resource_vcf? "--resource:dbsnp,known=true,training=false,truth=false,prior=2 ${dbsnp_resource_vcf} \\\n    ":"";

  my $faFile = get_option( $config, $section, "fasta_file", 1 );

  my $java_option = $classobj->get_java_option($config, $section, $memory);

  my $script = dirname(__FILE__) . "/fixLeftTrimDeletion.py";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

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

  print $pbs "
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
    $mills_option $axiomPoly_option ${dbsnp_option} -mode INDEL 

  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $task_name.MakeSitesOnlyVcf.failed
    rm $indels_recalibration
  fi
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

  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $task_name.SnpVariantRecalibrator.failed
    rm $snps_recalibration
  fi
fi
";

  my $rmlist = "$indels_recalibration ${indels_recalibration}.tbi \\
    $snps_recalibration ${snps_recalibration}.tbi";

  return($indels_recalibration, $indels_tranches, $snps_recalibration, $snps_tranches, $rmlist);
}

sub add_apply_vqsr_pbs {
  my ( $classobj, $config, $section, $pbs, $sample_name, $variant_filtered_vcf, $indels_recalibration, $indels_tranches, $snps_recalibration, $snps_tranches, $pass_file, $option, $thread, $memory ) = @_;

  my $java_option = $classobj->get_java_option($config, $section, $memory);
  
  my $indel_filter_level = get_option($config, $section, "indel_filter_level", 99.7);
  my $snp_filter_level = get_option($config, $section, "snp_filter_level", 99.7);

  my $reader_threads = min( 5, $thread );
  my $indel_recalibration_tmp_vcf = $sample_name . ".indels.recal.tmp.vcf.gz";
  my $recalibrated_vcf_filename = $sample_name . ".indels.snp.recal.vcf.gz";

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

  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sample_name.IndelApplyVQSR.failed
    rm $indel_recalibration_tmp_vcf
  fi
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

  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sample_name.SnpApplyVQSR.failed
    rm $recalibrated_vcf_filename
  fi
fi

if [[ -s $recalibrated_vcf_filename && ! -s $pass_file ]]; then
  echo SelectVariant=`date`
  gatk --java-options \"$java_option\" \\
    SelectVariants \\
    -O $pass_file \\
    -V $recalibrated_vcf_filename \\
    --exclude-filtered

  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sample_name.SelectVariant.failed
    rm $pass_file
  fi
fi

";

  my $rmlist = "$recalibrated_vcf_filename ${recalibrated_vcf_filename}.tbi \\
    $indel_recalibration_tmp_vcf ${indel_recalibration_tmp_vcf}.tbi";
  return ($rmlist);
}

sub add_left_trim_pbs {
  my ( $classobj, $config, $section, $pbs, $sample_name, $pass_file, $final_file ) = @_;

  my $script = dirname(__FILE__) . "/fixLeftTrimDeletion.py";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my $split_file = $sample_name . ".split.vcf";
  my $left_trim_file = $sample_name . ".norm.vcf";
  my $fix_file = $final_file;
  $fix_file =~ s/.gz//g;

print $pbs "
rm -f $split_file.failed $left_trim_file.failed

if [[ -s $pass_file && ! -s $left_trim_file ]]; then
  echo LeftAlignAndNorm=`date`
  bcftools norm -m- -o $split_file $pass_file 
  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $split_file.failed
    rm -f $split_file
  else
    bcftools norm -f $faFile -o $left_trim_file $split_file 
    status=\$?
    if [[ \$status -ne 0 ]]; then
      touch $left_trim_file.failed
      rm -f $left_trim_file
    else
      touch left_trim_file.succeed
      rm -f $split_file
    fi  
  fi  
fi

if [[ -s $left_trim_file && ! -s $final_file ]]; then
  echo noSpanDeletion=`date`
  python3 $script -i $left_trim_file -o $fix_file
  bgzip $fix_file
  tabix -p vcf $final_file
fi
";

  my $rmlist = "$split_file $left_trim_file";
  return ($rmlist);
}

1;
