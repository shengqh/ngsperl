#!/usr/bin/perl
package CQS::ValidFastqExtractor;

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
  $self->{_name}   = "ValidFastqExtractor";
  $self->{_suffix} = "_vfe";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $cqstools = get_cqstools( $config, $section, 1 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $sampleName ( sort keys %rawFiles ) {
    my $sampleFile = $rawFiles{$sampleName}->[0];
    my $finalFile  = $sampleName . "_valid.fastq.gz";

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

if [ -s $finalFile ]; then
  echo job has already been done. if you want to do again, delete ${resultDir}/${finalFile} and submit job again.
  exit 0;
fi

echo ValidFastqExtractor_start=`date`

mono $cqstools fastq_valid_extractor -i $sampleFile -o $finalFile 

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

  print "!!!shell file $shfile created, you can run this shell file to submit all bwa tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $extension = get_option( $config, $section, "extension" );
  my $merge_result = get_option( $config, $section, "merge_result", 0 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( sort keys %rawFiles ) {
    my $finalFile   = $resultDir . "/" . $sampleName . "_valid.fastq.gz";
    my @resultFiles = ();
    push( @resultFiles, $finalFile );
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
