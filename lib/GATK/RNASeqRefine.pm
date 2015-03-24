#!/usr/bin/perl
package GATK::RNASeqRefine;

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
  $self->{_name}   = "GATK::RNASeqRefine";
  $self->{_suffix} = "_rrf";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my @vcfFiles = @{ $config->{$section}{vcf_files} } or die "Define vcf_files in section $section first.";
  my $gatk_jar   = get_param_file( $config->{$section}{gatk_jar},   "gatk_jar",   1 );
  my $picard_jar = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1 );

  my $gatk_option = get_option( $config, $section, "gatk_option", "" );

  my $sorted = get_option( $config, $section, "sorted", 0 );

  my $knownvcf      = "";
  my $knownsitesvcf = "";

  foreach my $vcf (@vcfFiles) {
    $knownvcf      = $knownvcf . " -known $vcf";
    $knownsitesvcf = $knownsitesvcf . " -knownSites $vcf";
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles    = @{ $rawFiles{$sampleName} };
    my $sampleFile     = $sampleFiles[0];
    my $sampleFileName = basename($sampleFile);

    my $inputFile     = $sampleFile;
    my $presortedFile = "";
    my $sortCmd       = "";
    if ( !$sorted ) {
      my $presortedPrefix = $sampleName . ".sorted";
      $presortedFile = $presortedPrefix . ".bam";
      $sortCmd       = "samtools sort -@ $thread -m 4G $sampleFile $presortedPrefix";
      $inputFile     = $presortedFile;
    }

    my $rmdupFile = $sampleName . ".rmdup.bam";
    my $splitFile = $sampleName . ".rmdup.split.bam";
    my $grpFile   = $splitFile . ".grp";
    my $recalFile = $sampleName . ".rmdup.split.recal.bam";

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

    print SH "\$MYCMD ./$pbsName \n";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $curDir

echo GATKRNASeqRefine_start=`date` 

if [ -s $recalFile ]; then
  echo job has already been done. if you want to do again, delete $recalFile and submit job again.
  exit 0
fi

if [ ! -s $rmdupFile ]; then
  echo MarkDuplicates=`date` 
  $sortCmd
  java $option -jar $picard_jar MarkDuplicates I=$inputFile O=$rmdupFile ASSUME_SORTED=true REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${rmdupFile}.metrics
fi

if [ ! -s $splitFile ]; then
  echo SplitNCigarReads=`date` 
  java $option -jar $gatk_jar -T SplitNCigarReads -R $faFile -I $rmdupFile -o $splitFile -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -fixMisencodedQuals
fi

if [[ -s $splitFile && ! -s $grpFile ]]; then
  echo BaseRecalibrator=`date` 
  java $option -jar $gatk_jar -T BaseRecalibrator -rf BadCigar -R $faFile -I $splitFile $knownsitesvcf -o $grpFile
fi

if [[ -s $splitFile && -s $grpFile && ! -s $recalFile ]]; then
  echo PrintReads=`date`
  java $option -jar $gatk_jar -T PrintReads -rf BadCigar -R $faFile -I $splitFile -BQSR $grpFile -o $recalFile 
fi

if [[ -s $recalFile && ! -s ${recalFile}.bai ]]; then
  echo BamIndex=`date` 
  samtools index $recalFile
  samtools flagstat $recalFile > ${recalFile}.stat
  rm $presortedFile $rmdupFile ${sampleName}.rmdup.bai ${rmdupFile}.metrics $splitFile ${sampleName}.rmdup.split.bai $grpFile
fi
  
echo finished=`date`

exit 0;
";

    close OUT;

    print "$pbsFile created\n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all GATK refine tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $sortedFile  = $sampleName . ".rmdup.split.recal.bam";
    my @resultFiles = ();
    push( @resultFiles, "${resultDir}/${sampleName}/${sortedFile}" );
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
