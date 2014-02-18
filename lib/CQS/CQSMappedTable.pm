#!/usr/bin/perl
package CQS::CQSMappedTable;

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
  $self->{_name}   = "CQSMappedTable";
  $self->{_suffix} = "_mt";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqsFile = get_param_file( $config->{$section}{cqs_tools}, "cqs_tools", 1 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  my $shfile = $self->taskname( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "
cd $resultDir
";

  if ( defined $config->{$section}{groups} || defined $config->{$section}{groups_ref} ) {
    my $groups = get_raw_files( $config, $section, "groups" );
    for my $groupName ( sort keys %{$groups} ) {
      my @samples = @{ $groups->{$groupName} };

      my $filelist   = $self->getfile( $pbsDir,    "${task_name}_${groupName}", ".filelist", 0 );
      my $outputfile = $self->getfile( $resultDir, "${task_name}_${groupName}", ".count",    0 );
      my $outputname = basename($outputfile);

      open( FL, ">$filelist" ) or die "Cannot create $filelist";
      for my $sampleName ( sort @samples ) {
        my @countFiles = @{ $rawFiles{$sampleName} };
        my $countFile  = $countFiles[0];
        print FL $sampleName, "\t", $countFile, "\n";
      }
      close(FL);

      print SH "
cd $resultDir

mono-sgen $cqsFile mapped_table $option -o $outputname -l $filelist
";
    }
  }
  else {
    my $filelist   = $self->getfile( $pbsDir,    ${task_name}, ".filelist", 0 );
    my $outputfile = $self->getfile( $resultDir, ${task_name}, ".count",    0 );
    my $outputname = basename($outputfile);

    open( FL, ">$filelist" ) or die "Cannot create $filelist";
    for my $sampleName ( sort keys %rawFiles ) {
      my @countFiles = @{ $rawFiles{$sampleName} };
      my $countFile  = $countFiles[0];
      print FL $sampleName, "\t", $countFile, "\n";
    }
    close(FL);

    print SH "
cd $resultDir

mono-sgen $cqsFile mapped_table $option -o $outputname -l $filelist
";
  }
  close SH;
  print "!!!shell file $shfile created, you can run this shell file to run CQSMappedTable task.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $result = {};

  my @resultFiles = ();
  if ( defined $config->{$section}{groups} || defined $config->{$section}{groups_ref} ) {
    my $groups = get_raw_files( $config, $section, "groups" );
    for my $groupName ( sort keys %{$groups} ) {
      my $filelist   = $self->getfile( $pbsDir,    "${task_name}_${groupName}", ".filelist", 0 );
      my $outputfile = $self->getfile( $resultDir, "${task_name}_${groupName}", ".count",    0 );
      my $resultFile = "${resultDir}/$outputfile";

      push( @resultFiles, $resultFile );
      push( @resultFiles, $filelist );
    }
  }
  else {
    my $filelist   = $self->getfile( $pbsDir,    ${task_name}, ".filelist", 0 );
    my $outputfile = $self->getfile( $resultDir, ${task_name}, ".count",    0 );

    push( @resultFiles, $outputfile );
    push( @resultFiles, $filelist );
  }

  $result->{$task_name} = filter_array( \@resultFiles, $pattern );

  return $result;
}

1;
