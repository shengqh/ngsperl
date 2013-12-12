#!/usr/bin/perl
package CQS::ParalyzerClusterAnnotator;

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
  $self->{_name} = "CQS::ParalyzerClusterAnnotator";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqsFile   = get_param_file( $config->{$section}{cqs_tools},  "cqs_tools",  1 );
  my $gffFile   = get_param_file( $config->{$section}{gff_file},   "gff_file",   1 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}_an.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";

  for my $sampleName ( sort keys %rawFiles ) {
    my @bamFiles  = @{ $rawFiles{$sampleName} };
    my $bamFile   = $bamFiles[0];
    my $annFile = change_extension($bamFile, ".ann.csv");

    print SH "mono-sgen $cqsFile paralyzer_annotation $option -i $bamFile -g $gffFile -o $annFile
";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->name() . " tasks.\n";
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
