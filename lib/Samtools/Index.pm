#!/usr/bin/perl
package Samtools::Index;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::GroupTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "Samtools::Index";
  $self->{_suffix} = "_ix";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $isbamsorted = $config->{$section}{isbamsorted};
  if ( !defined($isbamsorted) ) {
    $isbamsorted = 0;
  }

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $pbsFile = $self->pbsfile($pbsDir, $sampleName);
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );
    
    print SH "\$MYCMD ./$pbsName \n";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

echo index=`date`
";

    my $bamFile = $sampleFiles[0];

    my $bamSortedFile;
    if ($isbamsorted) {
      $bamSortedFile = $bamFile;
    }
    else {
      ( $bamSortedFile, my $bamSorted ) = get_sorted_bam($bamFile);
      print OUT "if [ ! -s $bamSortedFile ]; then\n";
      print OUT "  echo samtools_sort=`date`\n";
      print OUT "  samtools sort $bamFile $bamSorted \n";
      print OUT "fi\n";
    }

    my $bamIndexFile = $bamSortedFile . ".bai";
    print OUT "if [ ! -s $bamIndexFile ]; then
  echo samtools_index=`date`
  samtools index $bamSortedFile 
fi

echo finished=`date`
";
    close OUT;

    print "$pbsFile created\n";
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $isbamsorted = $config->{$section}{isbamsorted};
  if ( !defined($isbamsorted) ) {
    $isbamsorted = 0;
  }

  my $result = {};
  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $bamFile = $sampleFiles[0];

    my $bamSortedFile;
    if ($isbamsorted) {
      $bamSortedFile = $bamFile;
    }
    else {
      ( $bamSortedFile, my $bamSorted ) = get_sorted_bam($bamFile);
    }

    my $bamIndexFile = $bamSortedFile . ".bai";

    my @resultFiles = ();
    push( @resultFiles, $bamIndexFile );
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }

  return $result;
}

1;
