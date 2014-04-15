#!/usr/bin/perl
package Bacteria::RockHopper;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::PairTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::PairTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "Bacteria::RockHopper";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $rawFiles = get_raw_files( $config, $section );
  my $groups   = get_raw_files( $config, $section, "groups" );
  my $pairs    = get_raw_files( $config, $section, "pairs" );
  my $rockhopper_jar = get_param_file( $config->{$section}{rockhopper_jar}, "rockhopper_jar", 1 );
  my $genome_dir  = get_option( $config, $section, "genome_dir",  1 );
  my $java_option = get_option( $config, $section, "java_option", "" );

  my %tpgroups = ();
  for my $groupName ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$groupName} };
    my @gfiles  = ();
    foreach my $sampleName ( sort @samples ) {
      my @fastqFiles = @{ $rawFiles->{$sampleName} };
      push( @gfiles, join( '%', @fastqFiles ) );
    }
    $tpgroups{$groupName} = join( ",", @gfiles );
  }

  my $mapfile = $self->getfile( $resultDir, $task_name, ".map" );
  open( MAP, ">$mapfile" ) or die "Cannot create $mapfile";
  print MAP "Fastq\tSample\tGroup\n";
  for my $groupName ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$groupName} };
    foreach my $sampleName ( sort @samples ) {
      my @fastqFiles = @{ $rawFiles->{$sampleName} };
      foreach my $fastq ( sort @fastqFiles ) {
        my $fname = basename($fastq);
        $fname = change_extension( $fname, "" );
        print MAP "$fname\t$sampleName\t$groupName\n";
      }
    }
  }
  close(MAP);

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct), "\n";

  for my $pairName ( sort keys %{$pairs} ) {
    my @groupNames = @{ $pairs->{$pairName} };
    if ( scalar(@groupNames) != 2 ) {
      die $pairName . " should include and only include two groups, currently is [" . join( ", ", @groupNames );
    }
    my @fastqs = ();
    push( @fastqs, $tpgroups{ $groupNames[1] } );
    push( @fastqs, $tpgroups{ $groupNames[0] } );
    my $fastqstrs = join( " ", @fastqs );

    my $pbsFile = $self->pbsfile( $pbsDir, $pairName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $pairName );

    my $curDir = create_directory_or_die( $resultDir . "/$pairName" );

    my $labels = $groupNames[1] . "," . $groupNames[0];

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

if [ -s ${curDir}/summary.txt ];then
  echo job has already been done. if you want to do again, delete ${curDir}/summary.txt and submit job again.
  exit 0;
fi

cd $curDir

java $java_option -cp $rockhopper_jar Rockhopper $option -g $genome_dir -L $labels -o . $fastqstrs

if [ -e intermediary ]; then
  rm -rf intermediary
fi

echo finished=`date`

exit 0
";

    close(OUT);

    print "$pbsFile created. \n";

    print SH "\$MYCMD ./$pbsName \n";
  }

  print SH "exit 0\n";
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all FastQC tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};

  my $rawFiles = get_raw_files( $config, $section );
  my $groups   = get_raw_files( $config, $section, "groups" );
  my $pairs    = get_raw_files( $config, $section, "pairs" );
  my $genome_name = get_option( $config, $section, "genome_name", 1 );

  for my $pairName ( sort keys %{$pairs} ) {
    my $curDir      = $resultDir . "/$pairName";
    my @resultFiles = ();
    push( @resultFiles, $curDir . "/" . $genome_name . "_operons.txt" );
    push( @resultFiles, $curDir . "/" . $genome_name . "_transcripts.txt" );
    $result->{$pairName} = filter_array( \@resultFiles, $pattern );
  }
  my $mapfile = $self->getfile( $resultDir, $task_name, ".map" );
  my @resultFiles = ();
  push( @resultFiles, $mapfile );
  $result->{$task_name} = filter_array( \@resultFiles, $pattern );

  return $result;
}

1;
