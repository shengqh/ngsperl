#!/usr/bin/perl
package Data::Wget;

use strict;
use warnings;
use File::Basename;
use List::Util;
use File::Slurp;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::Task;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "Wget";
  $self->{_suffix} = "_wget";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $batch = get_option( $config, $section, "batch", 20 );

  for my $sampleName ( sort keys %rawFiles ) {
    my @listFiles = @{ $rawFiles{$sampleName} };
    my $listFile  = $listFiles[0];
    my @urls      = read_file( $listFile, chomp => 1 );

    my $urlCount = scalar(@urls);

    for ( my $i = 0 ; $i < $urlCount ; $i += $batch ) {
      my $iend = $i + $batch;
      if ( $iend >= $urlCount ) {
        $iend = $urlCount;
      }

      my $name = "${sampleName}_${i}_${iend}";
      my $pbsFile = $self->pbsfile( $pbsDir, $name );
      my $pbsName = basename($pbsFile);
      my $log     = $self->logfile( $logDir, $name );

      my $log_desc = $cluster->get_log_desc($log);

      open( OUT, ">$pbsFile" ) or die $!;
      print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

echo Wget=`date`
";

      for ( my $j = $i ; $j < $iend ; $j++ ) {
        my $url = $urls[$j];
        my $filename = basename($url);
        
        print OUT "if [ ! -s $filename ]; then
  wget $url
fi

";
        
      }
       
      close(OUT);
    }
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result      = {};

  return $result;
}

1;
