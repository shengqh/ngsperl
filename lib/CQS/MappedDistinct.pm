#!/usr/bin/perl
package CQS::MappedDistinct;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::CombinedTask;

our @ISA = qw(CQS::CombinedTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "MappedDistinct";
  $self->{_suffix} = "_dt";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqsFile = get_param_file( $config->{$section}{cqs_tools}, "cqs_tools", 1 );

  my %firstFiles = %{ get_raw_files( $config, $section ) };
  my %secondFiles = %{ get_raw_files( $config, $section, "second" ) };
  my $firstSuffix  = $config->{$section}{first_suffix}  or die "define ${section}::first_suffix first";
  my $secondSuffix = $config->{$section}{second_suffix} or die "define ${section}::second_suffix first";

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "cd $resultDir
echo CQSMappedDistinct=`date`

";

  for my $sampleName ( sort keys %firstFiles ) {
    my $firstFile  = $firstFiles{$sampleName}->[0];
    my $secondFile = $secondFiles{$sampleName}->[0];
    my $firstoutput = $firstSuffix . $sampleName . ".distinct.count";
    my $secondoutput = $secondSuffix . $sampleName . ".distinct.count";
    
    print SH "mono-sgen $cqsFile mapped_distinct $option --inputfile1 $firstFile --outputfile1 $firstoutput --inputfile2 $secondFile --outputfile2 $secondoutput

";
  }
  print SH "
echo finished=`date`

exit 1
";
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all trna_count tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  my $firstSuffix  = $config->{$section}{first_suffix}  or die "define ${section}::first_suffix first";
  my $secondSuffix = $config->{$section}{second_suffix} or die "define ${section}::second_suffix first";

  my $result = {};
  for my $sampleName (sort keys %rawFiles ) {
    my $firstoutput = $resultDir . "/" . $firstSuffix . $sampleName . ".distinct.count";
    my $secondoutput = $resultDir . "/" . $secondSuffix . $sampleName . ".distinct.count";
    
    my @resultFiles = ();

    push( @resultFiles, $firstoutput );
    push( @resultFiles, $secondoutput );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
