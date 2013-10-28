#!/usr/bin/perl
package CQS::CQSDatatable;

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
  $self->{_name} = "CQSDatatable";
  bless $self, $class;
  return $self;
}

sub get_result{
  my ( $task_name, $option ) = @_;

  my $result;
  if($option =~ /-o\s+(\S+)/){
    $result = $1;
  }
  else{
    $result = $task_name . ".count";
    $option = $option . " -o " . $result;
  }
  return ($result, $option);
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqsFile = get_param_file( $config->{$section}{cqs_tools},     "cqs_tools",     1 );
  my $mapFile = get_param_file( $config->{$section}{name_map_file}, "name_map_file", 0 );

  my $mapoption = "";
  if (defined $mapFile) {
    $mapoption = "-m $mapFile";
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $filelist = $pbsDir . "/${task_name}.filelist";
  open( FL, ">$filelist" ) or die "Cannot create $filelist";
  for my $sampleName ( sort keys %rawFiles ) {
    my @bamFiles = @{ $rawFiles{$sampleName} };
    my $bamFile  = $bamFiles[0];
    print FL $sampleName, "\t", $bamFile, "\n";
  }
  close(FL);

  my ($result, $newoption) = get_result($task_name, $option);

  my $shfile = $pbsDir . "/${task_name}_tb.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "
cd $resultDir

mono-sgen $cqsFile data_table $newoption -l $filelist $mapoption
";

  close SH;

  print "!!!shell file $shfile created, you can run this shell file to run cqs_datatable task.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result     = {};
  my ($resultFile, $newoption) = get_result($task_name, $option);
  $resultFile = $resultDir . "/" . $resultFile;
  my $filelist   = $resultDir . "/${task_name}.filelist";

  my @resultFiles = ();
  push( @resultFiles, $resultFile );
  push( @resultFiles, $filelist );

  $result->{$task_name} = filter_array( \@resultFiles, $pattern );

  return $result;
}

1;
