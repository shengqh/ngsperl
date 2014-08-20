#!/usr/bin/perl
package Imputation::Impute2;

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
  $self->{_name}   = "Imputation::Impute2";
  $self->{_suffix} = "_imp";
  bless $self, $class;
  return $self;
}

sub containPosition {
  my ( $positions, $start, $end ) = @_;
  for my $position ( @{$positions} ) {
    if ( $position >= $start && $position <= $end ) {
      return 1;
    }
  }

  return 0;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  my $mergefile = $self->taskfile( $pbsDir, $task_name . "_merge" );
  open( MSH, ">$mergefile" ) or die "Cannot create $mergefile";

  my %rawFiles   = %{ get_raw_files( $config, $section ) };
  my %mapFiles   = %{ get_raw_files( $config, $section, "map_file" ) };
  my %haploFiles = %{ get_raw_files( $config, $section, "haplo_file" ) };
  my %pedFiles   = %{ get_raw_files( $config, $section, "ped_file" ) };

  my $maxChromosomeLength = get_option( $config, $section, "max_chromosome_length" );
  my $interval            = get_option( $config, $section, "interval" );

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $sample      = $sampleFiles[0];

    my @mapFiles = @{ $mapFiles{$sampleName} };
    my $map      = $mapFiles[0];

    my @haploFiles = @{ $haploFiles{$sampleName} };
    my $haploFile  = $haploFiles[0];
    my $legendFile = change_extension( $haploFile, ".legend" );

    my @peds   = @{ $pedFiles{$sampleName} };
    my $ped    = $peds[0];
    my $pedmap = change_extension( $ped, ".map" );

    my @positions;
    open( INFILE, "<", $pedmap ) or die("Couldn't open $pedmap for reading!\n");

    while (<INFILE>) {
      push @positions, ( split( /\s+/, $_ ) )[3];
    }

    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

    my $gen_file = "${sampleName}.gen";

    my $start = 1;
    while ( $start < $maxChromosomeLength ) {
      my $end = $start + $interval - 1;

      if ( !containPosition( \@positions, $start, $end ) ) {
        next;
      }

      my $cursample = $sampleName . "_" . $start;

      my $pbsFile = $self->pbsfile( $pbsDir, $cursample );
      my $pbsName = basename($pbsFile);
      my $log     = $self->logfile( $logDir, $cursample );

      my $tmpFile = "${cursample}.tmp";

      open( OUT, ">$pbsFile" ) or die $!;

      print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir 

if [ -s $gen_file ]; then
  echo job has already been done. if you want to do again, delete $curDir/$gen_file and submit job again.
  exit 0;
fi

if [ -s $tmpFile ]; then
  echo job has already been done. if you want to do again, delete $curDir/$tmpFile and submit job again.
  exit 0;
fi

echo impute2_start=`date` 

impute2 -known_haps_g $sample -m $map -int $start $end -h $haploFile -l $legendFile -o $tmpFile

echo finished=`date` 

";
      print MSH "cat $tmpFile >> $gen_file \n";
      $start = $end + 1;

      close(OUT);

      print SH "\$MYCMD ./$pbsName \n";
      print "$pbsFile created\n";
    }
  }
  print SH "exit 0\n";
  close(SH);

  close(MSH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my @resultFiles = ();

    push( @resultFiles, "${resultDir}/${sampleName}/${sampleName}.gen" );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
