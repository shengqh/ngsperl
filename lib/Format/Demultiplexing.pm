#!/usr/bin/perl
package Format::Demultiplexing;

use strict;
use warnings;
use File::Basename;
use Data::Table;
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
  $self->{_name}   = "Format::Demultiplexing";
  $self->{_suffix} = "_dem";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  my %mapFiles = %{ get_raw_files( $config, $section, "maps" ) };
  my $cqsFile = get_param_file( $config->{$section}{cqs_tools}, "cqs_tools", 1 );

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $fastqfile   = $sampleFiles[0];
    my $summaryfile = $sampleName . "_indecies.tsv";

    my @maps    = @{ $mapFiles{$sampleName} };
    my $mapfile = $maps[0];

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

echo demultiplexing_started=`date`

mono-sgen $cqsFile fastq_demultiplex $option -m $mapfile -i $fastqfile -o . -s $summaryfile

echo demultiplexing_finished=`date`

exit 0 
";
    close OUT;

    print "$pbsFile created \n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->name() . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %mapFiles = %{ get_raw_files( $config, $section, "maps" ) };

  my $result = {};
  for my $sampleName ( keys %mapFiles ) {
    my @maps    = @{ $mapFiles{$sampleName} };
    my $mapfile = $maps[0];

    my $table = Data::Table::fromTSV( $mapfile, 0 );
    if ( $table->nofcol() == 3 ) {
      foreach my $i ( 0 .. $table->lastRow ) {
        $result->{ $table->elm( $i, 3 ) } = $resultDir . "/" . $table->elm( $i, 2 );
      }
    }
    else {
      foreach my $i ( 0 .. $table->lastRow ) {
        my $name = change_extension_gzipped( $table->elm( $i, 2 ) );
        $result->{$name} = $resultDir . "/" . $table->elm( $i, 2 );
      }
    }
  }
  return $result;
}

1;
