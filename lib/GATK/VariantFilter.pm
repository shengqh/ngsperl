#!/usr/bin/perl
package GATK::VariantFilter;

use strict;
use warnings;
use File::Basename;
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $dbsnp = get_param_file( $config->{$section}{dbsnp_vcf}, "dbsnp_vcf", 1 );
  my $vqsrMode = get_option( $config, $section, "vqsr_mode" );
  my $hapmap;
  my $omni;
  my $g1000;
  my $mills;

  if ($vqsrMode) {
    $hapmap = get_param_file( $config->{$section}{hapmap_vcf}, "hapmap_vcf", 1 );
    $omni   = get_param_file( $config->{$section}{omni_vcf},   "omni_vcf",   0 );
    $g1000  = get_param_file( $config->{$section}{g1000_vcf},  "g1000_vcf",  0 );
    $mills  = get_param_file( $config->{$section}{mills_vcf},  "mills_vcf",  1 );
  }

  my $cqstools = get_cqstools( $config, $section, 1 );

  my $faFile   = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $gatk_jar = get_param_file( $config->{$section}{gatk_jar},   "gatk_jar",   1 );

  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option || $java_option eq "" ) {
    $java_option = "-Xmx${memory}";
  }

  my %gvcfFiles = %{ get_raw_files( $config, $section ) };

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $mergedFile = $task_name . ".raw.vcf";

  my $finalFile = $task_name . ".pass.vcf";

  my $log_desc = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $finalFile );

  print $pbs "
if [ ! -s $mergedFile ]; then
  echo GenotypeGVCFs=`date` 
  java $java_option -jar $gatk_jar -T GenotypeGVCFs $option -nt $thread -D $dbsnp -R $faFile \\
";

  for my $sample_name ( sort keys %gvcfFiles ) {
    my @sample_files = @{ $gvcfFiles{$sample_name} };
    my $gvcfFile     = $sample_files[0];
    print $pbs "    --variant $gvcfFile \\\n";
  }

  print $pbs "  -o $mergedFile
fi
";

  if ($vqsrMode) {    #VQSR mode
    my $snpCal        = $task_name . ".snp.recal";
    my $snpTranches   = $task_name . ".snp.tranche";
    my $snpPass       = $task_name . ".recal_snp_raw_indel.vcf";
    my $indelCal      = $task_name . ".indel.recal";
    my $indelTranches = $task_name . ".indel.tranche";
    my $indelPass     = $task_name . ".recal_snp_recal_indel.vcf";
    print $pbs "
if [[ -s $mergedFile && ! -s $snpCal ]]; then
  echo VariantRecalibratorSNP=`date` 
  java $java_option -jar $gatk_jar \\
    -T VariantRecalibrator -nt $thread \\
    -R $faFile \\
    -input $mergedFile \\
";

    if ($hapmap) {
      print $pbs "    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \\\n";
    }

    if ($omni) {
      print $pbs "    -resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \\\n";
    }

    if ($g1000) {
      print $pbs "    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $g1000 \\\n";
    }

    print $pbs "    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \\
    -an DP \\
    -an QD \\
    -an FS \\
    -an SOR \\
    -an MQ \\
    -an MQRankSum \\
    -an ReadPosRankSum \\
    -an InbreedingCoeff \\
    -mode SNP \\
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
    -recalFile $snpCal \\
    -tranchesFile $snpTranches
fi

if [[ -s $snpCal && ! -s $snpPass ]]; then
  echo ApplyRecalibrationSNP=`date` 
  java $java_option -jar $gatk_jar \\
    -T ApplyRecalibration -nt $thread \\
    -R $faFile \\
    -input $mergedFile \\
    -mode SNP \\
    --ts_filter_level 99.0 \\
    -recalFile $snpCal \\
    -tranchesFile $snpTranches \\
    -o $snpPass   
fi

if [[ -s $snpPass && ! -s $indelCal ]]; then
  echo VariantRecalibratorIndel=`date` 
  java $java_option -jar $gatk_jar \\
    -T VariantRecalibrator -nt $thread \\
    -R $faFile \\
    -input $snpPass \\
    -resource:mills,known=true,training=true,truth=true,prior=12.0 $mills \\
    -an QD \\
    -an DP \\
    -an FS \\
    -an SOR \\
    -an MQRankSum \\
    -an ReadPosRankSum \\
    -an InbreedingCoeff \\
    -mode INDEL \\
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
    --maxGaussians 4 \\
    -recalFile $indelCal \\
    -tranchesFile $indelTranches
fi

if [[ -s $snpPass && -s $indelCal && ! -s $indelPass ]]; then
  echo ApplyRecalibrationIndel=`date` 
  java $java_option -jar $gatk_jar \\
    -T ApplyRecalibration -nt $thread \\
    -R $faFile \\
    -input $snpPass \\
    -mode INDEL \\
    --ts_filter_level 99.0 \\
    -recalFile $indelCal \\
    -tranchesFile $indelTranches \\
    -o $indelPass
fi

if [[ -s $indelPass && ! -s $finalFile ]]; then
  cat $indelPass | awk '\$1 ~ \"#\" || \$7 == \"PASS\"' > $finalFile
  rm $snpCal ${snpCal}.idx $snpTranches $snpPass ${snpPass}.idx $indelCal ${indelCal}.idx $indelTranches $indelPass ${indelPass}.idx
fi


";
  }
  else {    #hard filter mode
    my $raw_snp_file   = $task_name . ".raw_snp.vcf";
    my $snpPass        = $task_name . ".snp_filtered.vcf";
    my $raw_indel_file = $task_name . ".raw_indel.vcf";
    my $indelPass      = $task_name . ".indel_filtered.vcf";
    my $filtered_file  = $task_name . ".filtered.vcf";

    my $snp_filter =
      get_option( $config, $section, "is_rna" )
      ? "-window 35 -cluster 3 --filterExpression \"FS > 30.0 || QD < 2.0\" --filterName \"snp_filter\" "
      : "--filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filterName \"snp_filter\" ";

    print $pbs "
if [[ -s $mergedFile && ! -s $snpPass ]]; then
  echo filter_snp=`date`
  java $java_option -Xmx${memory} -jar $gatk_jar -T SelectVariants -R $faFile -V $mergedFile -selectType SNP -o $raw_snp_file 
  java $java_option -Xmx${memory} -jar $gatk_jar -T VariantFiltration -R $faFile -V $raw_snp_file $snp_filter -o $snpPass 
fi
";

    my $indel_filter =
      ( get_option( $config, $section, "is_rna" ) ? "-window 35 -cluster 3" : "" ) . " --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" --filterName \"indel_filter\" ";

    print $pbs "
if [[ -s $mergedFile && ! -s $indelPass ]]; then
  echo filter_indel=`date`
  java $java_option -Xmx${memory} -jar $gatk_jar -T SelectVariants -R $faFile -V $mergedFile -selectType INDEL -o $raw_indel_file 
  java $java_option -Xmx${memory} -jar $gatk_jar -T VariantFiltration -R $faFile -V $raw_indel_file $indel_filter -o $indelPass 
fi

if [[ -s $snpPass && -s $indelPass && ! -s $finalFile ]]; then
  echo combine_snp=`date`
  java $java_option -Xmx${memory} -jar $gatk_jar -T CombineVariants -R $faFile --variant:snp $snpPass --variant:indel $indelPass -o $filtered_file -genotypeMergeOptions PRIORITIZE -priority snp,indel
  cat $filtered_file | awk '\$1 ~ \"#\" || \$7 == \"PASS\"' > $finalFile
  rm $raw_snp_file $raw_indel_file $snpPass $indelPass $filtered_file ${raw_snp_file}.idx ${raw_indel_file}.idx ${snpPass}.idx ${indelPass}.idx ${filtered_file}.idx
fi
"
  }
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $finalFile = $task_name . ".pass.vcf";
  my @result_files = ();
  push( @result_files, $result_dir . "/" . $finalFile );

  #  my @result_files = ();
  #  my $finalFile    = $task_name . ".median3.snp.pass.vcf";
  #  push( @result_files, $result_dir . "/" . $finalFile );
  #  $finalFile = $task_name . ".median3.indel.pass.vcf";
  #  push( @result_files, $result_dir . "/" . $finalFile );

  my $result = {};
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
