#!/usr/bin/perl
package CQS::Bam2Fastq;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "Bam2Fastq";
  $self->{_suffix} = "_b2q";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $ispaired = get_option( $config, $section, "ispaired", 0 );

  my $sort_before_convert = get_option( $config, $section, "sort_before_convert" );
  my $sort_thread         = get_option( $config, $section, "sort_thread" );
  my $sortoption = $sort_thread < 2 ? "" : "-@ $sort_thread";

  my $unmapped_only = get_option( $config, $section, "unmapped_only", 0 );
  my $unzipped      = get_option( $config, $section, "unzipped",      0 );
  if ($unzipped) {
    $option = $option . " -u";
  }

  print "sort_before_convert = $sort_before_convert\n";
  print "sort_thread = $sort_thread\n";

  my $cqstools = $config->{$section}{cqstools} or die "define ${section}::cqstools first";

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $bamfile     = $sampleFiles[0];

    my $pbsFile = $self->pbsfile($pbsDir, $sampleName);
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    open( OUT, ">$pbsFile" ) or die $!;

    my $finalFile = $ispaired ? $sampleName . ".1.fastq" : $sampleName . ".fastq";
    if ( !$unzipped ) {
      $finalFile = $finalFile . ".gz";
    }
    my $convertCmd = "";

    if ($sort_before_convert) {
      my $sourceFile = "${sampleName}.sortname.bam";
      my $sortCmd;
      if ($unmapped_only) {
        $sortCmd = "samtools view -b -f 4 $bamfile | samtools sort $option -n $sortoption - ${sampleName}.sortname";
      }
      else {
        $sortCmd = "samtools sort $option -n $sortoption $bamfile ${sampleName}.sortname";
      }

      $convertCmd = "if [ ! -s $sourceFile ]; then
    $sortCmd
  fi
  
  mono-sgen $cqstools bam2fastq $option -i $sourceFile -o $sampleName 
  
  if [ -s $finalFile ]; then
    rm $sourceFile
  fi";
    }
    else {
      if ($unmapped_only) {
        my $unmapped_bam = "${sampleName}.unmapped.bam";
        $convertCmd = "samtools view -b -f 4 $bamfile > $unmapped_bam
  mono-sgen $cqstools bam2fastq $option -i $unmapped_bam -o $sampleName
  rm $unmapped_bam ";
      }
      else {
        $convertCmd = "mono-sgen $cqstools bam2fastq $option -i $bamfile -o $sampleName ";
      }
    }

    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir

echo started=`date`

if [ ! -s $finalFile ]; then
  $convertCmd
fi

echo finished=`date`

exit 1 
";
    close OUT;

    print "$pbsFile created \n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all Bam2Fastq tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $ispaired = get_option( $config, $section, "ispaired", 0 );
  my $unzipped = get_option( $config, $section, "unzipped", 0 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {

    my @resultFiles = ();
    if ($ispaired) {
      if ($unzipped) {
        push( @resultFiles, $resultDir . "/" . $sampleName . ".1.fastq" );
        push( @resultFiles, $resultDir . "/" . $sampleName . ".2.fastq" );
      }
      else {
        push( @resultFiles, $resultDir . "/" . $sampleName . ".1.fastq.gz" );
        push( @resultFiles, $resultDir . "/" . $sampleName . ".2.fastq.gz" );
      }
    }
    else {
      if ($unzipped) {
        push( @resultFiles, $resultDir . "/" . $sampleName . ".fastq" );
      }
      else {
        push( @resultFiles, $resultDir . "/" . $sampleName . ".fastq.gz" );
      }
    }

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
