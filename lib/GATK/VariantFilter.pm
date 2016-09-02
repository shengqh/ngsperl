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
  my $vqsrMode = get_option( $config, $section, "vqsr_mode");
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

  my $snpCal      = $task_name . ".snp.recal";
  my $snpTranches = $task_name . ".snp.tranche";

  my $snpPass = $task_name . ".recal_snp_raw_indel.vcf";

  my $indelCal      = $task_name . ".indel.recal";
  my $indelTranches = $task_name . ".indel.tranche";

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

  if ( $vqsrMode ) {    #VQSR mode
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

if [[ -s $snpPass && -s $indelCal && ! -s $finalFile ]]; then
  echo ApplyRecalibrationIndel=`date` 
  java $java_option -jar $gatk_jar \\
    -T ApplyRecalibration -nt $thread \\
    -R $faFile \\
    -input $snpPass \\
    -mode INDEL \\
    --ts_filter_level 99.0 \\
    -recalFile $indelCal \\
    -tranchesFile $indelTranches \\
    -o $finalFile
  rm $snpCal ${snpCal}.idx $snpTranches $snpPass ${snpPass}.idx $indelCal ${indelCal}.idx $indelTranches 
fi
";
  }
  else {    #hard filter mode
    my $snp_filter =
      get_option( $config, $section, "is_rna" )
      ? "-window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\""
      : "\\
         --filterExpression \"QD < 2.0\" \\
         --filterName \"QD\" \\
         --filterExpression \"FS > 60.0\" \\
         --filterName \"FS\" \\
         --filterExpression \"MQ < 40\" \\
         --filterName \"MQ\" \\
         --filterExpression \"MQRankSum < -12.5\" \\
         --filterName \"MQRankSum\" \\
         --filterExpression \"ReadPosRankSum < -8.0\" \\
         --filterName \"ReadPosRankSum\"\\
         ";

    print $pbs "
if [[ -s $mergedFile && ! -s $snpPass ]]; then
  java $java_option -Xmx${memory} -jar $gatk_jar -T VariantFiltration -R $faFile -V $mergedFile $snp_filter -o $snpPass 
fi
";

    my $indelFilterOut = $task_name . ".snp_indel.filtered.vcf";
    my $indel_filter = ( get_option( $config, $section, "is_rna" ) ? "-window 35 -cluster 3" : "" ) . " \\
      --filterExpression \"QD < 2.0\" \\
      --filterName \"QD\" \\
      --filterExpression \"FS > 200.0\"\\
      --filterName \"FS\" \\
      --filterExpression \"ReadPosRankSum < -20.0\"\\
      --filterName \"ReadPosRankSum\"\\
      ";

    print $pbs "
if [[ -s $snpPass && ! -s $finalFile ]]; then
  java $java_option -Xmx${memory} -jar $gatk_jar -T VariantFiltration -R $faFile -V $snpPass $indel_filter -o $indelFilterOut 
  cat $indelFilterOut | awk '\$1 ~ \"#\" || \$7 == \"PASS\"' > $finalFile
  rm $snpPass ${snpPass}.idx $indelFilterOut ${indelFilterOut}.idx
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

  my $result = {};
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
