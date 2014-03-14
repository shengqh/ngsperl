#!/usr/bin/perl
package Homer::FindPeaks;

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
  $self->{_name}   = "Homer::FindPeaks";
  $self->{_suffix} = "_fp";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %tagDirectories = %{ get_raw_files( $config, $section ) };

  my $pairs = get_raw_files( $config, $section, "pairs" );

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  my $threadcount = get_pbs_thread( $config->{$section}{pbs} );

  for my $pairName ( sort keys %{$pairs} ) {
    my ( $ispaired, $gNames ) = get_pair_groups( $pairs, $pairName );
    
    my @groupNames = @{$gNames};

    my $controlTag = $tagDirectories{ $groupNames[0] }[0];
    my $sampleTag  = $tagDirectories{ $groupNames[1] }[0];

    my $pbsFile = $self->pbsfile( $pbsDir, $pairName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $pairName );

    my $finalFile = $pairName . ".tsv";

    open( OUT, ">$pbsFile" ) or die $!;

    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir 

echo homer_FindPeaks_start=`date` 

if [ -s $finalFile ];then
  echo job has already been done. if you want to do again, delete ${resultDir}/${finalFile} and submit job again.
  exit 0;
fi

findPeaks $sampleTag -i $controlTag -o $finalFile

echo homer_FindPeaks_finished=`date` 

exit 0

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

  my $pairs = get_raw_files( $config, $section, "pairs" );

  my $result = {};
  for my $pairName ( sort keys %{$pairs} ) {
    my @resultFiles = ();
    push( @resultFiles, "${resultDir}/${pairName}.tsv" );
    $result->{$pairName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
