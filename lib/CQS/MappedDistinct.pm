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
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "MappedDistinct";
  $self->{_suffix} = "_dt";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqsFile = get_cqstools( $config, $section, 1 );

  my %firstFiles = %{ get_raw_files( $config, $section ) };
  my %secondFiles = %{ get_raw_files( $config, $section, "second" ) };
  my $firstSuffix  = $config->{$section}{first_suffix}  or die "define ${section}::first_suffix first";
  my $secondSuffix = $config->{$section}{second_suffix} or die "define ${section}::second_suffix first";

  my $pbsFile = $self->pbsfile( $pbsDir, $task_name );
  my $log = $self->logfile( $logDir, $task_name );

  my $cluster = get_cluster( $config, $section );
  my $log_desc = $cluster->get_log_desc($log);

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
$log_desc

$path_file
cd $resultDir

echo CQSMappedDistinct=`date`

";

  for my $sampleName ( sort keys %firstFiles ) {
    my $firstFile    = $firstFiles{$sampleName}->[0];
    my $secondFile   = $secondFiles{$sampleName}->[0];
    my $firstoutput  = $firstSuffix . $sampleName . ".distinct.count";
    my $secondoutput = $secondSuffix . $sampleName . ".distinct.count";

    print OUT "mono-sgen $cqsFile mapped_distinct $option --inputfile1 $firstFile --outputfile1 $firstoutput --inputfile2 $secondFile --outputfile2 $secondoutput

";
  }
  print OUT "
echo finished=`date`

exit 0
";
  close(OUT);

  print "!!!pbs file $pbsFile created, you can run this shell file to run all MappedDistinct tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  my $firstSuffix  = $config->{$section}{first_suffix}  or die "define ${section}::first_suffix first";
  my $secondSuffix = $config->{$section}{second_suffix} or die "define ${section}::second_suffix first";

  my $result = {};
  for my $sampleName ( sort keys %rawFiles ) {
    my $firstoutput  = $resultDir . "/" . $firstSuffix . $sampleName . ".distinct.count";
    my $secondoutput = $resultDir . "/" . $secondSuffix . $sampleName . ".distinct.count";

    my @resultFiles = ();

    push( @resultFiles, $firstoutput );
    push( @resultFiles, $secondoutput );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
