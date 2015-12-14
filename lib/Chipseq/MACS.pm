#!/usr/bin/perl
package Chipseq::MACS;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::GroupTask;
use CQS::NGSCommon;
use Data::Dumper;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "Chipseq::MACS";
  $self->{_suffix} = "_macs";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %group_sample_map = %{ $self->get_group_sample_map( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $groupName ( sort keys %group_sample_map ) {
    my @sampleFiles = @{ $group_sample_map{$groupName} };
    my $sampleCount = scalar(@sampleFiles);

    if ( $sampleCount != 2 ) {
      die "SampleFile should be normal,tumor paired.";
    }

    my $curDir = create_directory_or_die( $resultDir . "/$groupName" );

    my $control = $sampleFiles[0][1];
    my $sample  = $sampleFiles[1][1];

    my $final  = "${groupName}_peaks.bed";

    my $pbsFile = $self->pbsfile( $pbsDir, $groupName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $groupName );

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

echo MACS_start=`date`

macs $option -t $sample -c $control -n $groupName

echo MACS_end=`date`

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
  for my $groupName ( sort keys %group_sample_map ) {
    my $curDir      = $resultDir . "/$groupName";
    my @resultFiles = ();
    push( @resultFiles, $curDir . "/${groupName}_peaks.bed" );
    push( @resultFiles, $curDir . "/${groupName}_MACS_wiggle/treat" );

    $result->{$groupName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
