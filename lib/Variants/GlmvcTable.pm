#!/usr/bin/perl
package Variants::GlmvcTable;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::UniqueTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "GlmvcTable";
  $self->{_suffix} = "_gt";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $glmvcfile = get_param_file( $config->{$section}{execute_file}, "execute_file", 1 );
  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $filelist = $self->getfile( $pbsDir, ${task_name}, ".filelist", 0 );
  open( FL, ">$filelist" ) or die "Cannot create $filelist";
  for my $sampleName ( sort keys %rawFiles ) {
    my @countFiles = @{ $rawFiles{$sampleName} };
    my $countFile  = $countFiles[0];
    print FL $sampleName, "\t", $countFile, "\n";
  }
  close(FL);

  my $pbsFile = $self->pbsfile( $pbsDir, $task_name );
  my $pbsName = basename($pbsFile);
  my $log     = $self->logfile( $logDir, $task_name );
  my $final   = $resultDir . "/" . $task_name . ".tsv";

  my $log_desc = $cluster->get_log_desc($log);

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
$log_desc

$path_file 

echo GlmvcTable=`date` 

cd $resultDir

if [ -s $final ]; then
  echo job has already been done. if you want to do again, delete $final and submit job again.
  exit 0
fi

mono $glmvcfile table -i $filelist -o $final
";

  close OUT;

  print "$pbsFile created \n";

  print "!!!pbs file $pbsFile created, you can run or submit this pbs file.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $final       = $resultDir . "/" . $task_name . ".tsv";
  my @resultFiles = ($final);

  my $result = { $task_name => filter_array( \@resultFiles, $pattern ) };
  return $result;
}
1;
