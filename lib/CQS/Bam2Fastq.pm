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
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $ispaired = get_option( $config, $section, "ispaired", 0 );

  my $sort_before_convert = get_option( $config, $section, "sort_before_convert" );
  my $sort_thread         = get_option( $config, $section, "sort_thread");
  my $sortoption = $sort_thread < 2 ? "" : "-@ $sort_thread";

  print "sort_before_convert = $sort_before_convert\n";
  print "sort_thread = $sort_thread\n";

  my $cqstools = $config->{$section}{cqstools} or die "define ${section}::cqstools first";

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}_b2q.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $bamfile     = $sampleFiles[0];

    my $pbsName = "${sampleName}_b2q.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${sampleName}_b2q.log";

    open( OUT, ">$pbsFile" ) or die $!;

    my $finalFile = $ispaired ? $sampleName . ".1.fastq.gz" : $sampleName . ".fastq.gz";
    my $convertCmd;

    if ($sort_before_convert) {
      my $sourceFile = "${sampleName}.sortname.bam";
      $convertCmd = "if [ ! -s $sourceFile ]; then
    samtools sort $option -n $sortoption $bamfile ${sampleName}.sortname
  fi
  
  mono-sgen $cqstools bam2fastq $option -i $sourceFile -o $sampleName 
  
  if [ -s $finalFile ]; then
    rm $sourceFile
  fi";
    }
    else {
      $convertCmd = "mono-sgen $cqstools bam2fastq $option -i $bamfile -o $sampleName ";
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

  my $ispaired = get_option( $config, $section, "ispaired" );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {

    my @resultFiles = ();
    if ($ispaired) {
      push( @resultFiles, $resultDir . "/" . $sampleName . ".1.fastq.gz" );
      push( @resultFiles, $resultDir . "/" . $sampleName . ".2.fastq.gz" );
    }
    else {
      push( @resultFiles, $resultDir . "/" . $sampleName . ".fastq.gz" );
    }

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
