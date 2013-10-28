#!/usr/bin/perl
package CQS::ReorderSam;

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
  $self->{_name} = "ReorderSam";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

    my $faFile         = get_param_file( $config->{$section}{fasta_file},     "fasta_file",     1 );
  my $jar            = get_param_file( $config->{$section}{jar},            "jar",            1 );

  my $shfile = $pbsDir . "/${task_name}_reorder.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $sampleFile  = $sampleFiles[0];

    my $pbsName = "${sampleName}_ro.pbs";
    my $pbsFile = $pbsDir . "/$pbsName";
    my $log     = $logDir . "/${sampleName}_ro.log";

    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );
    open( OUT, ">$pbsFile" ) or die $!;

    my $finalFile = "${sampleName}.reordered.bam";

    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir 

if [ -s $finalFile ]; then
  if [ ! -s ${finalFile}.bai ]; then
    samtools index $finalFile
    exit 1;
  fi
  echo job has already been done. if you want to do again, delete ${curDir}/${finalFile} and submit job again.
  exit 1;
fi

java -jar $jar I=${sampleFile} O=${finalFile} R=${faFile}

samtools index $finalFile

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
    my $finalFile = "${sampleName}.reordered.bam";
    my @resultFiles = ();
    push( @resultFiles, "${resultDir}/${sampleName}/${finalFile}" );
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
