#!/usr/bin/perl
package GATK::VariantFilterVQSR;

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
  $self->{_name}   = "VariantFilterVQSR";
  $self->{_suffix} = "_vf";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $dbsnp  = get_param_file( $config->{$section}{dbsnp_vcf},  "dbsnp_vcf",  1 );
  my $hapmap = get_param_file( $config->{$section}{hapmap_vcf}, "hapmap_vcf", 0 );
  my $omni   = get_param_file( $config->{$section}{omni_vcf},   "omni_vcf",   0 );
  my $g1000  = get_param_file( $config->{$section}{g1000_vcf},  "g1000_vcf",  0 );
  my $mills  = get_param_file( $config->{$section}{mills_vcf},  "mills_vcf",  0 );

  my $cqsFile = get_cqstools( $config, $section, 1 );

  my $faFile   = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $gatk_jar = get_param_file( $config->{$section}{gatk_jar},   "gatk_jar",   1 );

  my $min_median_depth = get_option( $config, $section, "min_median_depth", 3 );

  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option || $java_option eq "" ) {
    $java_option = "-Xmx${memory}";
  }

  my %gvcfFiles = %{ get_raw_files( $config, $section ) };

  my $pbsFile = $self->pbsfile( $pbsDir, $task_name );
  my $pbsName = basename($pbsFile);
  my $log     = $self->logfile( $logDir, $task_name );

  my $merged_file = $task_name . ".vcf";
  
  my $dpname = $task_name . ".median" . $min_median_depth;
  my $dpFilterOut = $dpname . ".vcf";

  my $snpOut      = $dpname . ".snp.vcf";
  my $snpCal      = $dpname . ".snp.cal";
  my $snpTranches = $dpname . ".snp.tranche";

  my $indelOut      = $dpname . ".indel.vcf";
  my $indelCal      = $dpname . ".indel.cal";
  my $indelTranches = $dpname . ".indel.tranche";

  my $snpPass   = $dpname . ".snp.pass.vcf";
  my $indelPass = $dpname . ".indel.pass.vcf";

  my $log_desc = $cluster->get_log_desc($log);

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

echo VariantFilterVQSR=`date` 

if [ ! -s $merged_file ]; then
  echo GenotypeGVCFs=`date` 
  java $java_option -jar $gatk_jar -T GenotypeGVCFs $option -nt $thread -D $dbsnp -R $faFile \\
";

  for my $sampleName ( sort keys %gvcfFiles ) {
    my @sampleFiles = @{ $gvcfFiles{$sampleName} };
    my $gvcfFile    = $sampleFiles[0];
    print OUT "    --variant $gvcfFile \\\n";
  }

  print OUT "  -o $merged_file
fi
";
  print OUT "
if [ ! -s $dpFilterOut ]; then
  echo VCF_MinimumMedianDepth_Filter=`date` 
  mono $cqsFile vcf_filter -i $merged_file -o $dpFilterOut -d $min_median_depth
fi 

if [[ -s $dpFilterOut && ! -s $snpOut ]]; then
  java $java_option -jar $gatk_jar -T SelectVariants -R $faFile -V $dpFilterOut -selectType SNP -o $snpOut 
  java $java_option -jar $gatk_jar -T SelectVariants -R $faFile -V $dpFilterOut -selectType INDEL -o $indelOut 
fi
";

  if ( $hapmap || $omni ) {
    print OUT "
if [[ -s $snpOut && ! -s $snpCal ]]; then
  echo VariantRecalibratorSNP=`date` 
  java $java_option -jar $gatk_jar \\
    -T VariantRecalibrator -nt $thread \\
    -R $faFile \\
    -input $snpOut \\
";

    if ($hapmap) {
      print OUT "    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \\\n";
    }

    if ($omni) {
      print OUT "    -resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \\\n";
    }

    if ($g1000) {
      print OUT "    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $g1000 \\\n";
    }

    print OUT "    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \\
    -an DP \\
    -an QD \\
    -an FS \\
    -an SOR \\
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
    -input $snpOut \\
    -mode SNP \\
    --ts_filter_level 99.0 \\
    -recalFile $snpCal \\
    -tranchesFile $snpTranches \\
    -o $snpPass   
fi
";
  }
  else {
    my $snpFilterOut = $task_name . ".snp.filtered.vcf";

    my $snp_filter =
      get_option( $config, $section, "is_rna" )
      ? "-window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\""
      : "--filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" -filterName \"snp_filter\"";

    print OUT "
if [[ -s $snpOut && ! -s $snpPass ]]; then
  java $java_option -Xmx${memory} -jar $gatk_jar -T VariantFiltration -R $faFile -V $snpOut $snp_filter -o $snpFilterOut 
  cat $snpFilterOut | awk '\$1 ~ \"#\" || \$7 == \"PASS\"' > $snpPass
fi
"
  }

  if ($mills) {
    print OUT "
if [[ -s $indelOut && ! -s $indelCal ]]; then
  echo VariantRecalibratorIndel=`date` 
  java $java_option -jar $gatk_jar \\
    -T VariantRecalibrator -nt $thread \\
    -R $faFile \\
    -input $indelOut \\
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

if [[ -s $indelOut && -s $indelCal && ! -s $indelPass ]]; then
  echo ApplyRecalibrationIndel=`date` 
  java $java_option -jar $gatk_jar \\
    -T ApplyRecalibration -nt $thread \\
    -R $faFile \\
    -input $indelOut \\
    -mode INDEL \\
    --ts_filter_level 99.0 \\
    -recalFile $indelCal \\
    -tranchesFile $indelTranches \\
    -o $indelPass 
fi
";
  }
  else {
    my $indelFilterOut = $task_name . ".indel.filtered.vcf";
    my $indel_filter =
      ( get_option( $config, $section, "is_rna" ) ? "-window 35 -cluster 3" : "" ) . " --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" -filterName \"indel_filter\"";

    print OUT "
if [[ -s $indelOut && ! -s $snpPass ]]; then
  java $java_option -Xmx${memory} -jar $gatk_jar -T VariantFiltration -R $faFile -V $indelOut $indel_filter -o $indelFilterOut 
  cat $indelFilterOut | awk '\$1 ~ \"#\" || \$7 == \"PASS\"' > $indelPass
fi
"
  }
  print OUT "
echo finished=`date`

exit 0
";
  close(OUT);

  print "$pbsFile created\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $min_median_depth = get_option( $config, $section, "min_median_depth", 3 );
  my $dpname = $task_name . ".median" . $min_median_depth;

  my $snpPass   = $dpname . ".snp.pass.vcf";
  my $indelPass = $dpname . ".indel.pass.vcf";

  my @resultFiles = ();
  push( @resultFiles, $resultDir . "/" . $snpPass );
  push( @resultFiles, $resultDir . "/" . $indelPass );

  my $result = {};
  $result->{$task_name} = filter_array( \@resultFiles, $pattern );

  return $result;
}

1;
