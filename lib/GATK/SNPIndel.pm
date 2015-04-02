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
  $self->{_suffix} = "_snp";
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

  my $rnafilter = get_option( $config, $section, "is_rna" ) ? "-window 35 -cluster 3 " : "";

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

  my %group_sample_map = %{ getGroupSampleMap( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  #print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $groupName ( sort keys %group_sample_map ) {
    my $curDir       = create_directory_or_die( $resultDir . "/$groupName" );
    my $listfilename = "${groupName}.list";
    my $listfile     = $curDir . "/$listfilename";
    open( LIST, ">$listfile" ) or die "Cannot create $listfile";
    my @sampleFiles = @{ $group_sample_map{$groupName} };
    foreach my $sampleFile (@sampleFiles) {
      print LIST $sampleFile . "\n";
    }
    close(LIST);

    my $snpOut       = $groupName . "_snp.vcf";
    my $snpStat      = $groupName . "_snp.stat";
    my $snpFilterOut = $groupName . "_snp_filtered.vcf";
    my $snpPass      = $groupName . "_snp_filtered.pass.vcf";

    my $pbsFile = $self->pbsfile( $pbsDir, $groupName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $groupName );

    print SH "\$MYCMD ./$pbsName \n";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $curDir

echo SNP=`date` 

if [ ! -s $snpOut ]; then
  java $java_option -jar $gatk_jar -T HaplotypeCaller $option -R $faFile -I $listfilename -stand_call_conf 20.0 -stand_emit_conf 20.0 -d dnsnp_file $compvcf --out $snpOut -dontUseSoftClippedBases -nct $thread
fi

if [ -s $snpOut ]; then
  java $java_option -Xmx${memory} -jar $gatk_jar -T VariantFiltration -R $faFile -V $snpOut $rnafilter -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -o $snpFilterOut 

  cat $snpFilterOut | awk '\$1 ~ \"#\" || \$7 == \"PASS\"' > $snpPass
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

  my %group_sample_map = %{ getGroupSampleMap( $config, $section ) };
  for my $groupName ( sort keys %group_sample_map ) {
    my $curDir      = $resultDir . "/$groupName";
    my $snpPass     = $groupName . "_snp_filtered.pass.vcf";
    my $indelPass   = $groupName . "_indel_filtered.pass.vcf";
    my @resultFiles = ();
    push( @resultFiles, "${curDir}/${snpPass}" );
    push( @resultFiles, "${curDir}/${indelPass}" );
    $result->{$groupName} = filter_array( \@resultFiles, $pattern );
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
      push( @resultFiles, $self->pbsfile( $pbsDir, $sampleName . "_snp" ) );
      push( @resultFiles, $self->pbsfile( $pbsDir, $sampleName . "_id" ) );
      $result->{$sampleName} = \@resultFiles;
    }
  }

  return $result;
}

1;
