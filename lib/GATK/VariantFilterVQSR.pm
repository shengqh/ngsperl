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

  my $dbsnp = get_param_file( $config->{$section}{dbsnp_vcf}, "dbsnp_vcf", 1 );
  my $hapmap = get_param_file( $config->{$section}{hapmap_vcf}, "hapmap_vcf", 1 );
  my $omni = get_param_file( $config->{$section}{omni_vcf}, "omni_vcf", 1 );
  my $g1000 = get_param_file( $config->{$section}{g1000_vcf}, "g1000_vcf", 1 );
  my $mills = get_param_file( $config->{$section}{mills_vcf}, "mills_vcf", 1 );

  my $faFile   = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $gatk_jar = get_param_file( $config->{$section}{gatk_jar},   "gatk_jar",   1 );

  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option || $java_option eq "" ) {
    $java_option = "-Xmx${memory}";
  }

  my %gvcfFiles = %{ get_raw_files( $config, $section ) };

  my $pbsFile = $self->pbsfile( $pbsDir, $task_name );
  my $pbsName = basename($pbsFile);
  my $log     = $self->logfile( $logDir, $task_name );

  my $merged_file = $task_name . ".vcf";
  my $recal_snp_file = $task_name . ".recal.snp.vcf";
  my $recal_snp_indel_file = $task_name . ".recalibrated_variants.vcf";
  my $recal_snp_indel_pass_file = $task_name . ".recalibrated_variants.pass.vcf";
  my $log_desc    = $cluster->get_log_desc($log);

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
  
if [[ -s $merged_file && ! -s recalibrate_SNP.recal ]]; then
  echo VariantRecalibratorSNP=`date` 
  java $java_option -jar $gatk_jar \\
    -T VariantRecalibrator -nt $thread \\
    -R $faFile \\
    -input $merged_file \\
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \\
    -resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \\
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $g1000 \\
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \\
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
    -recalFile recalibrate_SNP.recal \\
    -tranchesFile recalibrate_SNP.tranches \\
    -rscriptFile recalibrate_SNP_plots.R 
fi

if [[ -s recalibrate_SNP.recal && ! -s $recal_snp_file ]]; then
  echo ApplyRecalibrationSNP=`date` 
  java $java_option -jar $gatk_jar \\
    -T ApplyRecalibration -nt $thread \\
    -R $faFile \\
    -input $merged_file \\
    -mode SNP \\
    --ts_filter_level 99.0 \\
    -recalFile recalibrate_SNP.recal \\
    -tranchesFile recalibrate_SNP.tranches \\
    -o $recal_snp_file   
fi

if [[ -s $recal_snp_file && ! -s recalibrate_INDEL.recal ]]; then
  echo VariantRecalibratorIndel=`date` 
  java $java_option -jar $gatk_jar \\
    -T VariantRecalibrator -nt $thread \\
    -R $faFile \\
    -input $recal_snp_file \\
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
    -recalFile recalibrate_INDEL.recal \\
    -tranchesFile recalibrate_INDEL.tranches \\
    -rscriptFile recalibrate_INDEL_plots.R 
fi

if [[ -s $recal_snp_file && -s recalibrate_INDEL.recal && ! -s $recal_snp_indel_file ]]; then
  echo ApplyRecalibrationIndel=`date` 
  java $java_option -jar $gatk_jar \\
    -T ApplyRecalibration -nt $thread \\
    -R $faFile \\
    -input $recal_snp_file \\
    -mode INDEL \\
    --ts_filter_level 99.0 \\
    -recalFile recalibrate_INDEL.recal \\
    -tranchesFile recalibrate_INDEL.tranches \\
    -o $recal_snp_indel_file 
fi

if [[ -s $recal_snp_indel_file && ! -s $recal_snp_indel_pass_file ]]; then
  grep -e \"^#\" $recal_snp_indel_file > $recal_snp_indel_pass_file
  grep PASS $recal_snp_indel_file >> $recal_snp_indel_pass_file
fi

echo finished=`date`

exit 0
";

  close(OUT);

  print "$pbsFile created\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $recal_snp_indel_pass_file = $task_name . ".recalibrated_variants.pass.vcf";
  my $result = { $task_name => [ $resultDir . "/${recal_snp_indel_pass_file}" ] };

  return $result;
}

1;
