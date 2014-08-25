#!/usr/bin/perl
package CQS::CQSMirnaTable;

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
  $self->{_name}   = "CQSMirnaTable";
  $self->{_suffix} = "_mt";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqsFile = get_cqstools( $config, $section, 1 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $pbsFile = $self->pbsfile( $pbsDir, $task_name );
  my $pbsName = basename($pbsFile);
  my $log     = $self->logfile( $logDir, $task_name );

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir
";

  if ( defined $config->{$section}{groups} || defined $config->{$section}{groups_ref} ) {
    my $groups = get_raw_files( $config, $section, "groups" );
    for my $groupName ( sort keys %{$groups} ) {
      my $filelist   = $self->getfile( $pbsDir,    "${task_name}_${groupName}", ".filelist", 0 );
      my $outputfile = $self->getfile( $resultDir, "${task_name}_${groupName}", ".count",    0 );
      my $outputname = basename($outputfile);

      my @samples = @{ $groups->{$groupName} };
      open( FL, ">$filelist" ) or die "Cannot create $filelist";
      for my $sampleName ( sort @samples ) {
        my @countFiles = @{ $rawFiles{$sampleName} };
        my $countFile  = $countFiles[0];
        print FL $sampleName, "\t", $countFile, "\n";
      }
      close(FL);

      print OUT "
mono-sgen $cqsFile mirna_table $option -o $outputname -l $filelist
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

    print OUT "
mono-sgen $cqsFile mirna_table $option -o $outputname -l $filelist
";
  }
  close OUT;

  print "!!!shell file $pbsFile created.\n";
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
      my $outputfile = $self->getfile( $resultDir, "${task_name}_${groupName}", ".count",      0 );
      my $isomirfile = $self->getfile( $resultDir, "${task_name}_${groupName}", ".isomir.tsv", 0 );
      my $filelist   = $self->getfile( $pbsDir,    "${task_name}_${groupName}", ".filelist",   0 );
      push( @resultFiles, $outputfile );
      push( @resultFiles, $isomirfile );
      push( @resultFiles, $filelist );
    }
  }
  else {
    my $outputfile = $self->getfile( $resultDir, ${task_name}, ".count",      0 );
    my $isomirfile = $self->getfile( $resultDir, ${task_name}, ".isomir.tsv", 0 );
    my $filelist   = $self->getfile( $pbsDir,    ${task_name}, ".filelist",   0 );
    push( @resultFiles, $outputfile );
    push( @resultFiles, $isomirfile );
    push( @resultFiles, $filelist );
  }
  $result->{$task_name} = filter_array( \@resultFiles, $pattern );

  return $result;
}

1;
