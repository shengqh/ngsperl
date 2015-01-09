#!/usr/bin/perl
package Format::Bam2Fastq;

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
  $self->{_name}   = "Format::Bam2Fastq";
  $self->{_suffix} = "_b2q";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $ispaired = get_option( $config, $section, "ispaired", 0 );
  my $unzipped = get_option( $config, $section, "unzipped", 0 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $bamfile     = $sampleFiles[0];

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    my $finalFile = $ispaired ? $sampleName . "_1.fastq" : $sampleName . ".fastq";
    if ( !$unzipped ) {
      $finalFile = $finalFile . ".gz";
    }

    my $cluster = get_cluster( $config, $section );
    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

echo started=`date`

if [ ! -s $finalFile ]; then
  bam2fastq -o ${sampleName}#.fastq $bamfile
";
    if ( !$unzipped ) {
      if ($ispaired) {
        print OUT "  gzip ${sampleName}_1.fastq
  gzip ${sampleName}_2.fastq ";
      }
      else {
        print OUT "  gzip ${sampleName}.fastq";
      }
    }
    print OUT "
fi

echo finished=`date`

exit 0 
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
        push( @resultFiles, $resultDir . "/" . $sampleName . "_1.fastq" );
        push( @resultFiles, $resultDir . "/" . $sampleName . "_2.fastq" );
      }
      else {
        push( @resultFiles, $resultDir . "/" . $sampleName . "_1.fastq.gz" );
        push( @resultFiles, $resultDir . "/" . $sampleName . "_2.fastq.gz" );
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
