#!/usr/bin/perl
package CQS::AbstractBowtie;

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
  $self->{_name} = "AbstractBowtie";
  bless $self, $class;
  return $self;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $samonly = $config->{$section}{samonly};
  if ( !defined $samonly ) {
    $samonly = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( sort keys %rawFiles ) {
    my $samFile = $sampleName . ".sam";
    my $bamFile     = $sampleName . ".bam";

    my $curDir  = $resultDir . "/$sampleName";

    my $finalFile;
    if ($samonly) {
      $finalFile = $samFile;
    }
    else {
      $finalFile = $bamFile;
    }
    my @resultFiles = ();
    push( @resultFiles, $curDir . "/" . $finalFile );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }

  return $result;
}

1;
