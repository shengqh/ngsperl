#!/usr/bin/perl
package Chipseq::MACS2;

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
use Data::Dumper;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "Chipseq::MACS2";
  $self->{_suffix} = "_macs2";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

    my $samples  = join(" ", @sampleFiles);

    my $final  = "${sampleName}_peaks.bed";

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file 

cd $curDir

if [ -s $final ];then
  echo job has already been done. if you want to do again, delete ${curDir}/$final and submit job again.
  exit 0;
fi

echo MACS2_start=`date`

macs2 callpeak $option -t $samples -n $sampleName -B -q 0.01 --outdir . 

echo MACS2_end=`date`

exit 0
";

    close(OUT);

    print "$pbsFile created. \n";

    print SH "\$MYCMD ./$pbsName \n";
  }

  print SH "exit 0\n";
  close(SH);

  print "!!!shell file $shfile created, you can run this shell file to submit tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %group_sample_map = %{ $self->get_group_sample_map( $config, $section ) };

  my $result = {};
  for my $sampleName ( sort keys %group_sample_map ) {
    my $curDir      = $resultDir . "/$sampleName";
    my @resultFiles = ();
    push( @resultFiles, $curDir . "/${sampleName}_peaks.bed" );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
