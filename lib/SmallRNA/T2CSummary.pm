#!/usr/bin/perl
package SmallRNA::T2CSummary;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "SmallRNA::T2CSummary";
  $self->{_suffix} = "_t2c";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $cqsFile = get_cqstools( $config, $section, 1 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $filelist = $self->getfile( $pbsDir, ${task_name}, ".filelist", 0 );
  open( FL, ">$filelist" ) or die "Cannot create $filelist";
  for my $sampleName ( sort keys %rawFiles ) {
    my @countFiles = @{ $rawFiles{$sampleName} };
    my $countFile  = $countFiles[0];
    print FL $sampleName, "\t", $countFile, "\n";
  }
  close(FL);

  my $pbsFile  = $self->pbsfile( $pbsDir, $task_name );
  my $pbsName  = basename($pbsFile);
  my $log      = $self->logfile( $logDir, $task_name );
  my $log_desc = $cluster->get_log_desc($log);

  my $outputfile = $self->getfile( $resultDir, ${task_name}, "_t2c.tsv", 0 );
  my $outputname = basename($outputfile);

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

if [ -s $outputname ]; then
  echo job has already been done. if you want to do again, delete $outputfile and submit job again.
  exit 0;
fi

echo T2CSummary=`date`

mono $cqsFile smallrna_t2c_summary $option -o $outputname -l $filelist

echo T2CSummary_finished
";
  close OUT;

  print "!!!shell file $pbsFile created.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result = {};

  my @resultFiles = ();
  my $outputfile  = $self->getfile( $resultDir, ${task_name}, "_t2c.tsv", 0 );
  my $filelist    = $self->getfile( $pbsDir, ${task_name}, ".filelist", 0 );
  push( @resultFiles, $outputfile );
  push( @resultFiles, $filelist );
  $result->{$task_name} = filter_array( \@resultFiles, $pattern );

  return $result;
}

1;
