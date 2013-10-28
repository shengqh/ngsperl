#!/usr/bin/perl
package CQS::Sortbam;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "Sortbam";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $sort_by_query = get_option_value( $config->{$section}{sort_by_query}, 0 );

  my $shfile = $pbsDir . "/${task_name}_sort.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  my $threadcount = get_pbs_thread( $config->{$section}{pbs} );
  my %rawFiles = %{ get_raw_files( $config, $section ) };

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $sampleFile  = $sampleFiles[0];

    my $pbsName = "${sampleName}_sort.pbs";
    my $pbsFile = $pbsDir . "/$pbsName";
    my $log     = $logDir . "/${sampleName}_sort.log";

    open( OUT, ">$pbsFile" ) or die $!;

    my $sortcmd;
    my $finalPrefix;
    my $finalFile;
    if ($sort_by_query) {
      $finalPrefix = "${sampleName}.sortedname";
      $finalFile   = "${finalPrefix}.bam";
      $sortcmd     = "samtools sort $option -n -@ $threadcount $sampleFile $finalPrefix";
    }
    else {
      $finalPrefix = "${sampleName}.sorted";
      $finalFile   = "${finalPrefix}.bam";
      $sortcmd     = "samtools sort $option -@ $threadcount $sampleFile $finalPrefix
samtools index $finalFile";
    }

    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir 

if [ -s $finalFile ]; then
  echo job has already been done. if you want to do again, delete ${resultDir}/${finalFile} and submit job again.
  exit 1;
fi

$sortcmd

exit 1;
";
    close(OUT);

    print SH "\$MYCMD ./$pbsName \n";
    print "$pbsFile created\n";
  }
  print SH "exit 0\n";
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  my $sort_by_query = get_option_value( $config->{$section}{sort_by_query}, 0 );

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $finalFile;
    if ($sort_by_query) {
      $finalFile = "${sampleName}.sortedname.bam";
    }
    else {
      $finalFile = "${sampleName}.sorted.bam";
    }

    my @resultFiles = ();
    push( @resultFiles, "${resultDir}/${finalFile}" );
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
