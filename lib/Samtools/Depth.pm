#!/usr/bin/perl
package Samtools::Depth;

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
  $self->{_name} = "Samtools::Depth";
  $self->{_suffix} = "_dp";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %group_sample_map = %{ $self->get_group_sample_map( $config, $section ) };

  my $minimum_depth = $config->{$section}{minimum_depth};
  my $cqsFile;
  my $cqscommand = "";
  if(defined $minimum_depth){
    $cqsFile = get_cqstools( $config, $section, 1 );
    $cqscommand = " | mono $cqsFile depth_filter -d $minimum_depth"
  }

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $groupName ( sort keys %group_sample_map ) {
    my @sampleFiles = @{ $group_sample_map{$groupName} };
    my $sampleCount = scalar(@sampleFiles);
    my $samples = "";
    for (my $index = 0; $index < $sampleCount; $index ++) {
      $samples = $samples . " " . $sampleFiles[$index][1];
    }

    my $depth = "${groupName}.depth";

    my $pbsFile = $self->pbsfile($pbsDir, $groupName);
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $groupName );

    print SH "\$MYCMD ./$pbsName \n";

    my $cluster = get_cluster( $config, $section );
    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file 

echo depth=`date` 

cd $resultDir

if [ ! -s $depth ]; then
    samtools depth $option $samples $cqscommand > $depth
fi

echo finished=`date`
";
    close OUT;

    print "$pbsFile created \n";
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

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $groupName ( keys %{$groups} ) {
    my @resultFiles = ();
    my $curDir      = $resultDir . "/$groupName";
    my $depth = "${groupName}.depth";
    push( @resultFiles, "$curDir/$depth" );
    $result->{$groupName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
