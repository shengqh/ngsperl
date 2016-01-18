#!/usr/bin/perl
package SmallRNA::SamExtract;

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
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_fm";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $cqsFile = get_cqstools( $config, $section, 1 );
  my %countFiles = %{ get_raw_files( $config, $section ) };
  my %bamFiles = %{ get_raw_files( $config, $section, "bam_files" ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  my $moduleName = __PACKAGE__;

  for my $sampleName ( sort keys %countFiles ) {
    my @sampleCounts = @{ $countFiles{$sampleName} };
    my $sampleCount  = $sampleCounts[0];

    my @sampleBams = @{ $bamFiles{$sampleName} };
    my $sampleBam  = $sampleBams[0];

    my $finalSam  = $sampleName . ".mapped.sam";
    my $finalFile = $sampleName . ".mapped.bam";

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
  echo job has already been done. if you want to do again, delete $finalFile and submit job again.
  exit 0
fi

echo ${moduleName}_start=`date` 

mono $cqsFile sam_extract $option --bam $sampleBam --count $sampleCount -o $finalSam

if [ -s $finalSam ]; then
  samtools view -S -b $finalSam > $finalFile
  if [ -s $finalFile ]; then
    samtools index $finalFile
    rm $finalSam
  fi
fi
  
echo ${moduleName}_end=`date`

exit 0 
";

    close OUT;

    print "$pbsFile created \n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all ${moduleName} tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );
  my %rawFiles = %{ get_raw_files( $config, $section ) };
  my $result = {};
  for my $sampleName ( sort keys %rawFiles ) {
    my @resultFiles = ();
    push( @resultFiles, $resultDir . "/" . $sampleName . ".mapped.bam" );
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
