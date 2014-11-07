#!/usr/bin/perl
package Imputation::Shapeit;

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
  $self->{_name}   = "Imputation::Shapeit";
  $self->{_suffix} = "_si";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  my %mapFiles = %{ get_raw_files( $config, $section, "genetic_map_file" ) };

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $sample      = $sampleFiles[0];
    
    my @mapFiles = @{ $mapFiles{$sampleName} };
    my $map = $mapFiles[0];

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );
    my $source_option;

    if ( $sample =~ /bed$/ ) {
      my $fam = change_extension( $sample, ".fam" );
      my $bim = change_extension( $sample, ".bim" );
      $source_option = "-B $sample $bim $fam";
    }
    elsif ( $sample =~ /ped$/ ){
      my $ped_map = change_extension( $sample, ".map" );
      $source_option = "-P $sample $ped_map";
    }
    else{
      my $gen_sample = change_extension( $sample, ".sample" );
      $source_option = "-G $sample $gen_sample";
    }

    my $phase_file  = "${sampleName}.haps";
    my $sample_file = "${sampleName}.sample";

    open( OUT, ">$pbsFile" ) or die $!;

    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir 

if [ -s $phase_file ]; then
  echo job has already been done. if you want to do again, delete ${curDir}/${phase_file} and submit job again.
  exit 0;
fi

echo shapeit_start=`date` 

shapeit $option $source_option -M $map -O $phase_file $sample_file

echo finished=`date` 

";
    close(OUT);

    print SH "\$MYCMD ./$pbsName \n";
    print "$pbsFile created\n";
  }
  print SH "exit 0\n";
  close(SH);

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

    push( @resultFiles, "${resultDir}/${sampleName}/${sampleName}.haps" );
    push( @resultFiles, "${resultDir}/${sampleName}/${sampleName}.sample" );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
