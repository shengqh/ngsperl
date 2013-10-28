#!/usr/bin/perl
package CQS::GATKSNPIndel;

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
  $self->{_name} = "GATKSNPIndel";
  bless $self, $class;
  return $self;
}

sub getGroupSampleMap{
  my ($config, $section ) = @_;
  
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

  return (\%group_sample_map)
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $faFile   = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $gatk_jar = get_param_file( $config->{$section}{gatk_jar},   "gatk_jar",   1 );

  my @vcfFiles = @{ $config->{$section}{vcf_files} };
  my $knownvcf = "";
  foreach my $vcf (@vcfFiles) {
    if ( $knownvcf eq "" ) {
      $knownvcf = "-D $vcf";
    }
    else {
      $knownvcf = $knownvcf . " -comp $vcf";
    }
  }

  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option ) {
    $java_option = "";
  }

  my $filter_snp_option = $config->{$section}{filter_snp_option};
  if ( !defined $filter_snp_option ) {
    $filter_snp_option =
"--filterExpression \"QD<2.0\" --filterName \"QD\" --filterExpression \"MQ<40.0\" --filterName \"MQ\" --filterExpression \"FS >60.0\" --filterName \"FS\" --filterExpression \"HaplotypeScore >13.0\" --filterName \"HaplotypeScore\" --filterExpression \"MQRankSum<-12.5\" --filterName \"MQRankSum\" --filterExpression \"ReadPosRankSum<-8.0\" --filterName \"ReadPosRankSum\"";
  }

  my $filter_indel_option = $config->{$section}{filter_indel_option};
  if ( !defined $filter_indel_option ) {
    $filter_indel_option =
"--filterExpression \"QD<2.0\" --filterName \"QD\" --filterExpression \"ReadPosRankSum<-20.0\" --filterName \"ReadPosRankSum\" --filterExpression \"InbreedingCoeff < -0.8\" --filterName \"InbreedingCoeff\" --filterExpression \"FS > 200.0\" --filterName \"FS\"";
  }

  my %group_sample_map = %{getGroupSampleMap($config, $section)};
  my $shfile = $pbsDir . "/${task_name}_snp.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

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

    my $pbsName = "${groupName}_snp.pbs";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${groupName}_snp.log";

    my $pbsFile = "${pbsDir}/$pbsName";
    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "
$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir

echo SNP=`date` 

if [ ! -s $snpOut ]; then
  java $java_option -jar $gatk_jar -T UnifiedGenotyper $option -R $faFile -I $listfilename $knownvcf --out $snpOut -metrics $snpStat -glm SNP
fi

java $java_option -jar $gatk_jar -T VariantFiltration $filter_snp_option -R $faFile -o $snpFilterOut --variant $snpOut 

cat $snpFilterOut | awk '\$1 ~ \"#\" || \$7 == \"PASS\"' > $snpPass 

echo finished=`date`
";
    close OUT;
    print "$pbsFile created\n";

    my $indelOut         = $groupName . "_indel.vcf";
    my $indelStat        = $groupName . "_indel.stat";
    my $indelFilteredOut = $groupName . "_indel_filtered.vcf";
    my $indelPass        = $groupName . "_indel_filtered.pass.vcf";

    $pbsName = "${groupName}_id.pbs";

    print SH "\$MYCMD ./$pbsName \n";

    $log = "${logDir}/${groupName}_id.log";

    $pbsFile = "${pbsDir}/$pbsName";
    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "
$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir

echo InDel=`date` 

if [ ! -s $indelOut ]; then
  java -jar $java_option $gatk_jar $option -T UnifiedGenotyper -R $faFile -I $listfilename $knownvcf --out $indelOut -metrics $indelStat -glm INDEL
fi

java $java_option -jar $gatk_jar -T VariantFiltration $filter_indel_option -R $faFile -o $indelFilteredOut --variant $indelOut

cat $indelFilteredOut | awk '\$1 ~ \"#\" || \$7 == \"PASS\"' > $indelPass 

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

  my %group_sample_map = %{getGroupSampleMap($config, $section)};
  for my $groupName ( sort keys %group_sample_map ) {
    my $curDir      = $resultDir . "/$groupName";
    my $snpPass     = $groupName . "_snp_filtered.pass.vcf";
    my $indelPass   = $groupName . "_indel_filtered.pass.vcf";
    my @resultFiles = ();
    push( @resultFiles, "${curDir}/${snpPass}" );
    push( @resultFiles, "${curDir}/${indelPass}" );
    $result->{$groupName} = filter_array(\@resultFiles, $pattern);
  }
  return $result;
}

1;
