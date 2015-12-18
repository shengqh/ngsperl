#!/usr/bin/perl
package Chipseq::MACS2Callpeak;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::GroupTask;
use CQS::NGSCommon;
use CQS::StringUtils;
use Data::Dumper;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "Chipseq::MACS2Callpeak";
  $self->{_suffix} = "_mc";
  bless $self, $class;
  return $self;
}

sub get_current_raw_files {
  my ( $self, $config, $section ) = @_;
  my $rawFiles;
  if ( has_raw_files( $config, $section, "groups" ) ) {
    $rawFiles = $self->get_group_samplefile_map( $config, $section );
  }
  else {
    $rawFiles = get_raw_files( $config, $section );
  }
  return $rawFiles;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %rawFiles = %{ $self->get_current_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $curDir      = create_directory_or_die( $resultDir . "/$sampleName" );

    my $samples = join( " ", @sampleFiles );

    my $final = "${sampleName}_peaks.bed";

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file 

cd $curDir

echo MACS2_start=`date`

if [ ! -s ${sampleName}_peaks.narrowPeak ]; then
  macs2 callpeak $option -t $samples -n $sampleName
fi

if [ ! -s 


echo MACS2_end=`date`

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

  print "!!!shell file $shfile created, you can run this shell file to submit tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ $self->get_current_raw_files( $config, $section ) };
  my $result = {};
  for my $sampleName ( sort keys %rawFiles ) {
    my $curDir      = $resultDir . "/$sampleName";
    my @resultFiles = ();
    push( @resultFiles, $curDir . "/${sampleName}_treat_pileup.bdg" );
    push( @resultFiles, $curDir . "/${sampleName}_control_lambda.bdg" );
    push( @resultFiles, $curDir . "/${sampleName}_peaks.narrowPeak" );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
