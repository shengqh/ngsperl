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
use CQS::UniqueTask;
use CQS::StringUtils;

our @ISA = qw(CQS::AbstractBuildSummary);

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

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $proteomicstools = get_param_file( $config->{$section}{proteomicstools}, "proteomicstools", 1 );

  my ( $datasets, $lines, $dataset ) = get_datasets( $config, $section );
  my %datasets = %{$datasets};
  my @lines    = @{$lines};
  my @dataset  = @{$dataset};

  my $currentParamFile = $resultDir . "/" . $task_name . ".param";
  open( OUT, ">$currentParamFile" ) or die $!;

  for ( my $index = 0 ; $index < scalar(@lines) ; $index++ ) {
    print OUT $lines[$index] . "\n";
    if ( $lines[$index] =~ "<Datasets>" ) {
      last;
    }
  }

  #print @dataset;

  for my $sampleName ( sort keys %datasets ) {
    print OUT "    <Dataset>\n";
    print OUT "      <Name>$sampleName</Name>\n";
    foreach my $dsline (@dataset){
      print OUT $dsline . "\n";
    }
    print OUT "      <PathNames>\n";
    my @sampleFiles = @{ $datasets{$sampleName} };
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
  
  my $resultFile = $resultDir . "/" . $task_name . ".noredundant";

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

echo buildsummary=`date`

if [ ! -s $resultFile ]; then
  mono $proteomicstools buildsummary -i $currentParamFile 
fi

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
