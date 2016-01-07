#!/usr/bin/perl
package Visualization::Coverage;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::GroupTask;
use CQS::StringUtils;

our @ISA = qw(CQS::GroupTask);

my $directory;

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "Visualization::Coverage";
  $self->{_suffix} = "_vc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $bedFiles = get_raw_files( $config, $section );
  my $groups   = get_raw_files( $config, $section, "groups" );
  my $bamFiles = get_raw_files( $config, $section, "bam_files" );

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $name ( sort keys %{$bedFiles} ) {
    my $curDir = create_directory_or_die( $resultDir . "/$name" );

    my @curBedFiles = @{ $bedFiles->{$name} };
    my @curBamNames = @{ $groups->{$name} };
    my @curBamFiles = ();
    for my $bamName (@curBamNames) {
      push( @curBamFiles, $bamFiles->{$bamName}->[0] );
    }

    my $curBamFileStr = join( ' ', @curBamFiles );

    my $bamCount = scalar(@curBamNames);

    my $configFileName = "${name}.filelist";
    my $configFile     = $curDir . "/${configFileName}";
    open( CON, ">$configFile" ) or die "Cannot create $configFile";
    print CON "Name\tFile\n";
    my $cutindecies = "1,2";
    my $curcutindex = 1;
    for ( my $index = 0 ; $index < $bamCount ; $index++ ) {
      print CON $curBamNames[$index], "\t", $curBamFiles[$index], "\n";
      $curcutindex += 3;
      $cutindecies = $cutindecies . "," . $curcutindex;
    }
    close CON;

    my $pbsFile = $self->pbsfile( $pbsDir, $name );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $name );

    my $log_desc = $cluster->get_log_desc($log);

    #my $final = "${comparisonName}_c3.0_common.bed";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file 

echo coverage=`date` 

cd $curDir

";

    for my $bedFile (@curBedFiles) {
      open( IIN, $bedFile ) or die $!;
      while (<IIN>) {
        s/\r|\n//g;
        my ( $chr, $start, $end ) = split "\t";
        print OUT "
  samtools mpileup -r ${chr}:${start}-${end} $curBamFileStr | cut -f${$cutindecies} > ${chr}_${start}_${end}.mpileup
";
      }
      close IIN;
    }

    print OUT "

echo finished=`date`

";
    close OUT;

    print "$pbsFile created \n";

    print SH "\$MYCMD ./$pbsName \n";
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

  my $comparisons = get_raw_files( $config, $section, "groups" );
  my $result = {};
  for my $comparisonName ( sort keys %{$comparisons} ) {
    my @resultFiles = ();
    push( @resultFiles, $resultDir . "/${comparisonName}_c3.0_common.bed" );
    push( @resultFiles, $resultDir . "/${comparisonName}_cond1.bed" );
    push( @resultFiles, $resultDir . "/${comparisonName}_cond2.bed" );
    $result->{$comparisonName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
