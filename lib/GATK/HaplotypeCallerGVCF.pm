#!/usr/bin/perl
package GATK::HaplotypeCallerGVCF;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "GATK::HaplotypeCallerGVCF";
  $self->{_suffix} = "_hc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $faFile   = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $gatk_jar = get_param_file( $config->{$section}{gatk_jar},   "gatk_jar",   1 );

  my $dbsnp   = get_param_file( $config->{$section}{dbsnp_vcf},      "dbsnp_vcf",      1 );
  my $compvcf = get_param_file( $config->{$section}{comparison_vcf}, "comparison_vcf", 0 );

  if ( defined $compvcf ) {
    $compvcf = "-comp " . $compvcf;
  }
  else {
    $compvcf = "";
  }

  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option ) {
    $java_option = "-Xmx${memory}";
  }

  my %bamFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  #print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $sampleName ( sort keys %bamFiles ) {
    my @sampleFiles = @{ $bamFiles{$sampleName} };
    my $bamFile = $sampleFiles[0];

    my $curDir       = create_directory_or_die( $resultDir . "/$sampleName" );

    my $snvOut       = $sampleName . "_snv.gvcf";
    my $snvStat      = $sampleName . "_snv.stat";

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $curDir

echo HaplotypeCaller=`date`

if [ ! -s $snvOut ]; then
  java $java_option -jar $gatk_jar -T HaplotypeCaller $option -R $faFile -I $bamFile -D $dbsnp $compvcf -nct $thread --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 --out $snvOut
fi

echo finished=`date`
";
    close OUT;
    print "$pbsFile created\n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all GATK SnpInDel tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );
  my $result = {};

  my %bamFiles = %{ get_raw_files( $config, $section ) };
  for my $sampleName ( sort keys %bamFiles ) {
    my $curDir      = $resultDir . "/$sampleName";
    my $snvOut      = $sampleName . "_snv.gvcf";
    my @resultFiles = ();
    push( @resultFiles, "${curDir}/${snvOut}" );
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
