#!/usr/bin/perl
package SmallRNA::FastqTrna;

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
  $self->{_name}   = "FastqTrna";
  $self->{_suffix} = "_ft";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $cqsFile = get_cqstools( $config, $section, 1 );
  my %rawFiles = %{ get_raw_files( $config, $section ) };
  my $extension = get_option( $config, $section, "extension" );

  my %seqCountFiles = ();
  if ( has_raw_files( $config, $section, "seqcount" ) ) {
    %seqCountFiles = %{ get_raw_files( $config, $section, "seqcount" ) };
  }

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $sampleFile  = $sampleFiles[0];
    my $finalFile1  = $sampleName . "_40less" . $extension;
    my $finalFile2  = $sampleName . "_40plus" . $extension;
    my $summaryFile = $sampleName . $extension . ".summary";

    my $seqcountFile = "";
    if ( defined $seqCountFiles{$sampleName} ) {
      my @seqcounts = @{ $seqCountFiles{$sampleName} };
      my $seqcount  = $seqcounts[0];
      $seqcountFile = " -c $seqcount";
    }

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

if [ -s $finalFile1 ]; then
  echo job has already been done. if you want to do again, delete $finalFile1 and $finalFile2 and submit job again.
  exit 0
fi

echo FastqTrna=`date` 

mono-sgen $cqsFile fastq_trna $option -i $sampleFile -l $finalFile1 -p $finalFile2 -s $summaryFile $seqcountFile

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

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  my $extension = get_option( $config, $section, "extension" );

  my %seqCountFiles = ();
  if ( defined $config->{$section}{"seqcount"} || defined $config->{$section}{"seqcount_ref"} ) {
    %seqCountFiles = %{ get_raw_files( $config, $section, "seqcount" ) };
  }

  my $result = {};
  for my $sampleName ( sort keys %rawFiles ) {
    my $finalFile1  = $resultDir . "/" . $sampleName . ".40less" . $extension;
    my $finalFile2  = $resultDir . "/" . $sampleName . ".40plus" . $extension;
    my $summaryFile = $resultDir . "/" . $sampleName . $extension . ".summary";

    my @resultFiles = ();
    push( @resultFiles, $finalFile1 );
    if ( defined $seqCountFiles{$sampleName} ) {
      push( @resultFiles, $finalFile1 . ".dupcount" );
    }
    push( @resultFiles, $finalFile2 );
    if ( defined $seqCountFiles{$sampleName} ) {
      push( @resultFiles, $finalFile2 . ".dupcount" );
    }
    push( @resultFiles, $summaryFile );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
