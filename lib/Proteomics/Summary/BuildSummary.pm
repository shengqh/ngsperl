#!/usr/bin/perl
package Proteomics::Summary::BuildSummary;

use strict;
use warnings;
use File::Basename;
use File::Slurp;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "Proteomics::Summary::BuildSummary";
  $self->{_suffix} = "_bs";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $proteomicstools = get_param_file( $config->{$section}{proteomicstools}, "proteomicstools", 1 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $parameterFile = get_param_file( $config->{$section}{parameter_file}, "parameter_file", 1 );

  my @lines = read_file( $parameterFile, chomp => 1 );

  my $currentParamFile = $pbsDir . "/" . $task_name . ".param";
  open( OUT, ">$currentParamFile" ) or die $!;
  my @dataset   = ();
  my $indataset = 0;
  for ( my $index = 0 ; $index < scalar(@lines) ; $index++ ) {
    if ( $lines[$index] =~ "<Dataset>" ) {
      $indataset = 1;
    }
    elsif ( $lines[$index] =~ "</Dataset>" ) {
      $indataset = 0;
    }
    elsif ( $indataset && $lines[$index] !~ "PathName" ) {
      push( @dataset, $lines[$index] );
    }
  }

  print @dataset;

  my $pbsFile = $self->pbsfile( $pbsDir, $task_name );
  my $pbsName = basename($pbsFile);
  my $log     = $self->logfile( $logDir, $task_name );

  my $log_desc = $cluster->get_log_desc($log);

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
$log_desc

$path_file
";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
  }

  print OUT "
echo finished=`date`

exit 0 
";
  close OUT;

  print "$pbsFile created \n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  my $rank2 = ( $option =~ /--rank2/ ) && ( $option =~ /Comet/ );
  for my $sampleName ( keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my @resultFiles = ();

    for my $sampleFile (@sampleFiles) {
      my $resultFile = $rank2 ? $sampleFile . ".rank2.peptides" : $sampleFile . ".peptides";
      push( @resultFiles, "${resultFile}" );
    }
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
