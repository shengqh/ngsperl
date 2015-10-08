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

  my @dataset   = ();
  my $indataset = 0;
  for ( my $index = 0 ; $index < scalar(@lines) ; $index++ ) {
    if ( $lines[$index] =~ "<Dataset>" ) {
      $indataset = 1;
    }
    elsif ( $lines[$index] =~ "</Dataset>" ) {
      $indataset = 0;
    }
    elsif ( !$indataset ) {
      next;
    }
    elsif ( $lines[$index] =~ "PathName" ) {
      next;
    }
    elsif ( $lines[$index] =~ "<Name>" ) {
      next;
    }
    else {
      push( @dataset, $lines[$index] );
    }
  }

  my $currentParamFile = $pbsDir . "/" . $task_name . ".param";
  open( OUT, ">$currentParamFile" ) or die $!;

  for ( my $index = 0 ; $index < scalar(@lines) ; $index++ ) {
    print OUT $lines[$index] . "\n";
    if ( $lines[$index] =~ "<Datasets>" ) {
      last;
    }
  }

  #print @dataset;

  for my $sampleName ( sort keys %rawFiles ) {
    print OUT "    <Dataset>\n";
    print OUT "      <Name>$sampleName</Name>\n";
    foreach my $dsline (@dataset){
      print OUT $dsline . "\n";
    }
    print OUT "      <PathNames>\n";
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    for my $sampleFile (@sampleFiles) {
      print OUT "        <PathName>$sampleFile</PathName>\n";
    }
    print OUT "      </PathNames>\n";
    print OUT "    </Dataset>\n";
  }

  for ( my $index = 0 ; $index < scalar(@lines) ; $index++ ) {
    if ( $lines[$index] =~ "</Datasets>" ) {
      for ( my $nextindex = $index ; $nextindex < scalar(@lines) ; $nextindex++ ) {
        print OUT $lines[$nextindex] . "\n";
      }
      last;
    }
  }

  close(OUT);

  my $pbsFile = $self->pbsfile( $pbsDir, $task_name );
  my $pbsName = basename($pbsFile);
  my $log     = $self->logfile( $logDir, $task_name );

  my $log_desc = $cluster->get_log_desc($log);

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
$log_desc

$path_file

echo buildsummary=`date`

mono $proteomicstools buildsummary -i $currentParamFile 

echo finished=`date`

exit 0 
";
  close OUT;

  print "$pbsFile created \n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result           = {};
  my $currentParamFile = $pbsDir . "/" . $task_name . ".param";
  my @resultFiles      = ();
  push( @resultFiles, $currentParamFile );
  $result->{$task_name} = filter_array( \@resultFiles, $pattern );
  return $result;
}

1;
