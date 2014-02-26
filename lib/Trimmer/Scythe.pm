#!/usr/bin/perl
package Trimmer::Scythe;

use strict;
use warnings;
use File::Basename;
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
  $self->{_name}   = "Trimmer::Scythe";
  $self->{_suffix} = "_scy";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $selfname = $self->{_name};

  my $faFile = get_param_file( $config->{$section}{adapter_file}, "adapter_file", 1 );
  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir

echo scythe_start=`date` 

";

    if ( scalar(@sampleFiles) == 2 ) {
      my $sample1 = $sampleFiles[0];
      my $sample2 = $sampleFiles[1];

      my $trim1 = change_extension_gzipped(basename($sample1), "_scythe.fastq");
      my $trim2 = change_extension_gzipped(basename($sample2), "_scythe.fastq");

      my $finalFile1 = $trim1 . ".gz";
      my $finalFile2 = $trim2 . ".gz";

      print OUT "
if [ ! -s $finalFile1 ]; then
  if [ ! -s $trim1 ]; then
    scythe $option -a $faFile -o $trim1 $sample1
  fi

  if [ ! -s $trim2 ]; then
    scythe $option -a $faFile -o $trim2 $sample2
  fi

  gzip $trim1
  gzip $trim2
fi
";
    }
    else {
      my $sample1 = $sampleFiles[0];

      my $trim1 = change_extension_gzipped(basename($sample1), "_scythe.fastq");

      my $finalFile1 = $trim1 . ".gz";

      print OUT "
if [ ! -s $finalFile1 ]; then
  if [ ! -s $trim1 ]; then
    scythe $option -a $faFile -o $trim1 $sample1
  fi

  gzip $trim1
fi
";
    }
    print OUT "
echo finished=`date`

exit 1 
";
    close OUT;

    print "$pbsFile created \n";
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my @resultFiles = ();

    if ( scalar(@sampleFiles) == 2 ) {
      my $sample1 = $sampleFiles[0];
      my $sample2 = $sampleFiles[1];

      my $trim1 = change_extension_gzipped(basename($sample1), "_scythe.fastq");
      my $trim2 = change_extension_gzipped(basename($sample2), "_scythe.fastq");

      my $finalFile1 = $trim1 . ".gz";
      my $finalFile2 = $trim2 . ".gz";

      push( @resultFiles, "${resultDir}/${finalFile1}" );
      push( @resultFiles, "${resultDir}/${finalFile2}" );
    }
    else {
      my $sample1 = $sampleFiles[0];

      my $trim1 = change_extension_gzipped(basename($sample1), "_scythe.fastq");

      my $finalFile1 = $trim1 . ".gz";

      push( @resultFiles, "${resultDir}/${finalFile1}" );
    }
    
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
