#!/usr/bin/perl
package GATK::GenotypeGVCFs;

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
  $self->{_name}   = "GenotypeGVCFs";
  $self->{_suffix} = "_gg";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $dbsnp = get_param_file( $config->{$section}{dbsnp_vcf}, "dbsnp_vcf", 0 );
  my $dbsnpparam = "";
  if ( defined $dbsnp ) {
    $dbsnpparam = "-D $dbsnp";
  }

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

  my $merged_file = $task_name . ".gvcf";
  my $log_desc    = $cluster->get_log_desc($log);

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

echo GenotypeGVCFs=`date` 

java $java_option -jar $gatk_jar -T GenotypeGVCFs $option -R $faFile \\
";

  for my $sampleName ( sort keys %gvcfFiles ) {
    my @sampleFiles = @{ $gvcfFiles{$sampleName} };
    my $gvcfFile    = $sampleFiles[0];
    print OUT "  --variant $gvcfFile \\\n";
  }

  print OUT "  -o $merged_file

echo finished=`date`

exit 0
";

  close(OUT);

  print "$pbsFile created\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $merged_file = $task_name . ".gvcf";
  my $result = { $task_name => [ $resultDir . "/${merged_file}" ] };

  return $result;
}

1;
