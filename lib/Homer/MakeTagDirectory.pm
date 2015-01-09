#!/usr/bin/perl
package Homer::MakeTagDirectory;

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
  $self->{_name}   = "Homer::MakeTagDirectory";
  $self->{_suffix} = "_mtd";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $rawFiles = get_raw_files( $config, $section );
  my $groups = get_raw_files( $config, $section, "groups" );
  
  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  my $threadcount = get_pbs_thread( $config->{$section}{pbs} );

  for my $groupName ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$groupName} };
    my @gfiles  = ();
    foreach my $sampleName ( sort @samples ) {
      my @bamFiles = @{ $rawFiles->{$sampleName} };
      foreach my $bam (@bamFiles){
        push( @gfiles, $bam );
      }
    }
    my $bams = join( " ", @gfiles );

    my $pbsFile = $self->pbsfile( $pbsDir, $groupName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $groupName );

    my $curDir      = create_directory_or_die( $resultDir . "/$groupName" );

    my $cluster = get_cluster( $config, $section );
    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $curDir 

echo homer_MakeTagDirectory_start=`date` 

makeTagDirectory $curDir $bams

echo homer_MakeTagDirectory_finished=`date` 

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

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $groupName ( sort keys %{$groups} ) {
    my @resultFiles = ();
    push( @resultFiles, "${resultDir}/${groupName}" );
    $result->{$groupName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
