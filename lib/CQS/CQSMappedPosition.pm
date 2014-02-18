#!/usr/bin/perl
package CQS::CQSMappedPosition;

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
  $self->{_name} = "CQS::CQSMappedPosition";
  $self->{_suffix} = "_mp";
  bless $self, $class;
  return $self;
}

sub get_result {
  my ( $task_name, $option ) = @_;

  my $result;
  if ( $option =~ /-o\s+(\S+)/ ) {
    $result = $1;
  }
  else {
    $result = $task_name . ".position";
    $option = $option . " -o " . $result;
  }
  return ( $result, $option );
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqsFile = get_param_file( $config->{$section}{cqs_tools}, "cqs_tools", 1 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  
  my $filelist = $pbsDir . "/" . $self->getname($task_name, ".filelist");
  open( FL, ">$filelist" ) or die "Cannot create $filelist";
  for my $sampleName ( sort keys %rawFiles ) {
    my @bamFiles = @{ $rawFiles{$sampleName} };
    my $bamFile  = $bamFiles[0];
    print FL $sampleName, "\t", $bamFile, "\n";
  }
  close(FL);

  my ( $result, $newoption ) = get_result( $task_name, $option );

  my $pbsFile = $self->pbsfile( $pbsDir, $task_name );
  my $pbsName = basename($pbsFile);
  my $log     = $self->logfile( $logDir, $task_name );

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir

mono-sgen $cqsFile mapped_position $newoption -l $filelist
";

  close OUT;

  print "!!!shell file $pbsFile created.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result = {};
  my ( $resultFile, $newoption ) = get_result( $task_name, $option );
  $resultFile = $resultDir . "/" . $resultFile;

  my $filelist = $pbsDir . "/" . $self->getname($task_name, ".filelist");

  my @resultFiles = ();
  push( @resultFiles, $resultFile );
  push( @resultFiles, $filelist );

  $result->{$task_name} = filter_array( \@resultFiles, $pattern );

  return $result;
}

1;
