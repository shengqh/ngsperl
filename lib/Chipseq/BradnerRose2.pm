#!/usr/bin/perl
package Chipseq::BradnerRose2;

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
  $self->{_name}   = "Chipseq::BradnerRose2";
  $self->{_suffix} = "_br";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %group_sample_map = %{ $self->get_group_sample_map( $config, $section ) };
  my $pipeline_dir = get_directory( $config, $section, "pipeline_dir", 1 );
  my $binding_site_file = get_param_file( $config, $section, "binding_site_file", 1 );

  if ( $option == "" ) {
    $option = "-s 12500 -t 2500";
  }

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

    my $normal = $sampleFiles[0][1];
    my $tumor  = $sampleFiles[1][1];

    my $cpRawFile  = "${groupName}.copynumber";
    my $cpCallFile = "${groupName}.call";
    my $cpSegFile  = "${groupName}.segment";

    my $pbsFile = $self->pbsfile( $pbsDir, $groupName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $groupName );

    print SH "\$MYCMD ./$pbsName \n";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file 

if [ -s gene_exp.diff ];then
  echo job has already been done. if you want to do again, delete ${curDir}/gene_exp.diff and submit job again.
  exit 0;
fi

echo BradnerRose2_start=`date`

cd $pipeline_dir
python ROSE2_main.py -i $binding_site_file -r $tumor -c $normal -o $curDir $option

echo end=`date`

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
    $result->{$groupName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
