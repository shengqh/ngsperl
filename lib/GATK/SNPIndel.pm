#!/usr/bin/perl
package GATK::SNPIndel;

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
  $self->{_name}   = "GATK::SNPIndel";
  $self->{_suffix} = "_snv";
  bless $self, $class;
  return $self;
}

sub getGroupSampleMap {
  my ( $config, $section ) = @_;

  my $rawFiles = get_raw_files( $config, $section );
  my %group_sample_map = ();
  if ( defined $config->{$section}{groups} || defined $config->{$section}{groups_ref} ) {
    my $groups = get_raw_files( $config, $section, "groups" );
    for my $groupName ( sort keys %{$groups} ) {
      my @samples = @{ $groups->{$groupName} };
      my @gfiles  = ();
      my $index   = 0;
      foreach my $sampleName (@samples) {
        my @bamFiles = @{ $rawFiles->{$sampleName} };
        push( @gfiles, $bamFiles[0] );
      }
      $group_sample_map{$groupName} = \@gfiles;
    }
  }
  else {
    %group_sample_map = %{$rawFiles};
  }

  return ( \%group_sample_map );
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $faFile   = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $gatk_jar = get_param_file( $config->{$section}{gatk_jar},   "gatk_jar",   1 );

  my $dbsnp   = get_param_file( $config->{$section}{dbsnp_vcf},      "dbsnp_vcf",      1 );
  my $compvcf = get_param_file( $config->{$section}{comparison_vcf}, "comparison_vcf", 0 );

  my $call_option = get_option( $config, $section, "is_rna" ) ? "-stand_emit_conf 10 -stand_call_conf 30":"-stand_emit_conf 20 -stand_call_conf 20"; 
  my $snp_filter = get_option( $config, $section, "is_rna" ) ? "-window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\"" : "--filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" -filterName \"snp_filter\"";  
  my $indel_filter = get_option( $config, $section, "is_rna" ) ? "-window 35 -cluster 3 --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\"" : "--filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\"";  

  if ( defined $compvcf ) {
    $compvcf = "-comp " . $compvcf;
  }
  else {
    $compvcf = "";
  }

  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option ) {
    $java_option = "";
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

    my $snvOut       = $sampleName . "_snv.vcf";
    my $snvStat      = $sampleName . "_snv.stat";

    my $snpOut       = $sampleName . "_snp.vcf";
    my $snpStat      = $sampleName . "_snp.stat";
    my $snpFilterOut = $sampleName . "_snp_filtered.vcf";
    my $snpPass      = $sampleName . "_snp_filtered.pass.vcf";

    my $indelOut       = $sampleName . "_indel.vcf";
    my $indelStat      = $sampleName . "_indel.stat";
    my $indelFilterOut = $sampleName . "_indel_filtered.vcf";
    my $indelPass      = $sampleName . "_indel_filtered.pass.vcf";

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

echo SNP=`date` 

if [[ -s $snpOut && ! -s $indelOut ]]; then
  mv $snpOut $snvOut
  if [ -s $snpStat ]; then
    mv $snpStat $snvStat
  fi 
fi

if [ ! -s $snvOut ]; then
  java $java_option -jar $gatk_jar -T HaplotypeCaller $option --genotyping_mode DISCOVERY -dontUseSoftClippedBases $call_option -R $faFile -I $bamFile -D $dbsnp $compvcf --out $snvOut -nct $thread 
fi

if [ -s $snvOut ]; then
  java $java_option -Xmx${memory} -jar $gatk_jar -T SelectVariants -R $faFile -V $snvOut -selectType SNP -o $snpOut 
  java $java_option -Xmx${memory} -jar $gatk_jar -T VariantFiltration -R $faFile -V $snpOut $snp_filter -o $snpFilterOut 
  cat $snpFilterOut | awk '\$1 ~ \"#\" || \$7 == \"PASS\"' > $snpPass

  java $java_option -Xmx${memory} -jar $gatk_jar -T SelectVariants -R $faFile -V $snvOut -selectType INDEL -o $indelOut 
  java $java_option -Xmx${memory} -jar $gatk_jar -T VariantFiltration -R $faFile -V $indelOut $indel_filter -o $indelFilterOut 
  cat $indelFilterOut | awk '\$1 ~ \"#\" || \$7 == \"PASS\"' > $indelPass
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
    my $snpPass      = $sampleName . "_snp_filtered.pass.vcf";
    my $indelPass      = $sampleName . "_indel_filtered.pass.vcf";
    my @resultFiles = ();
    push( @resultFiles, "${curDir}/${snpPass}" );
    push( @resultFiles, "${curDir}/${indelPass}" );
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

sub pbsfiles {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result = {};
  if ( $self->{_pbskey} eq "" ) {
    $result->{$task_name} = $self->pbsfile( $pbsDir, $task_name );
  }
  else {
    my %group_sample_map = %{ getGroupSampleMap( $config, $section ) };

    for my $sampleName ( sort keys %group_sample_map ) {
      my @resultFiles = ();
      push( @resultFiles, $self->pbsfile( $pbsDir, $sampleName ) );
      $result->{$sampleName} = \@resultFiles;
    }
  }

  return $result;
}

1;
