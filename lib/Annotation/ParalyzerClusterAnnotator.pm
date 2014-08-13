#!/usr/bin/perl
package Annotation::ParalyzerClusterAnnotator;

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
  $self->{_name} = "CQS::ParalyzerClusterAnnotator";
  $self->{_suffix} = "_an";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqsFile   = get_cqstools( $config, $section, 1 );
  my $corFiles   = $config->{$section}{coordinate_files} or die "define coordinate_files (array) in section $section first!";

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->pbsfile( $pbsDir, $task_name );
  my $log = $self->logfile( $logDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file
cd $resultDir

";

  for my $sampleName ( sort keys %rawFiles ) {
    my @bamFiles  = @{ $rawFiles{$sampleName} };
    my $bamFile   = $bamFiles[0];
    my $annFile = change_extension($bamFile, ".ann.csv");
    
    my $cfiles = merge_string(',', @{$corFiles});

    print SH "mono-sgen $cqsFile paralyzer_annotation $option -i $bamFile -c $cfiles -o $annFile
";
  }
  close(SH);

  print "!!!pbs file $shfile created.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $fasta_format = $config->{$section}{fasta_format};
  if ( !defined $fasta_format ) {
    $fasta_format = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my @bamFiles  = @{ $rawFiles{$sampleName} };
    my $bamFile   = $bamFiles[0];
    my $annFile = change_extension($bamFile, ".ann.csv");

    my @resultFiles = ();
    push( @resultFiles, $annFile );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
