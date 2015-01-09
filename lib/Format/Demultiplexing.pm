#!/usr/bin/perl
package Format::Demultiplexing;

use strict;
use warnings;
use File::Basename;
use Data::Table;
use String::Util qw(trim);
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
  my $cqsFile = get_cqstools( $config, $section, 1 );

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

    my $cluster = get_cluster( $config, $section );
    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

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
    if ( $table->nofCol() == 3 ) {
      foreach my $i ( 0 .. $table->lastRow ) {
        my $name     = trim( $table->elm( $i, 2 ) );
        my $filename = trim( $table->elm( $i, 1 ) );
        $result->{$name} = [ $resultDir . "/" . $filename ];
      }
    }
    else {
      foreach my $i ( 0 .. $table->lastRow ) {
        my $filename = trim( $table->elm( $i, 1 ) );
        my $name = change_extension_gzipped($filename);
        $result->{$name} = [ $resultDir . "/" . $filename ];
      }
    }
  }
  return $result;
}

1;
