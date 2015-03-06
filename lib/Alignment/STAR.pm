#!/usr/bin/perl
package Alignment::STAR;

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
  $self->{_name}   = "STAR";
  $self->{_suffix} = "_star";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  if ( $option !~ /outSAMprimaryFlag/ ) {
    $option = $option . "  --outSAMprimaryFlag AllBestScore";
  }

  my $sort_by_coordinate = get_option( $config, $section, "sort_by_coordinate", 0 );
  my $output_format = $sort_by_coordinate ? "--outSAMtype BAM SortedByCoordinate" : "--outSAMtype BAM Unsorted";

  my $genome_dir = $config->{$section}{genome_dir} or die "define ${section}::genome_dir first";
  my %fqFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  my $threadcount = get_pbs_thread( $config->{$section}{pbs} );

  for my $sampleName ( sort keys %fqFiles ) {
    my @sampleFiles = @{ $fqFiles{$sampleName} };

    my $uncompress = ( $sampleFiles[0] =~ /.gz$/ ) ? " --readFilesCommand zcat" : "";

    my $samples = join( " ", @sampleFiles );

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );
    my $rgline = "–-outSAMattrRGline ID:$sampleName SM:$sampleName LB:$sampleName";

    my $final = $sort_by_coordinate ? $sampleName . "_Aligned.sortedByCoord.out.bam" : $sampleName . "_Aligned.out.bam";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $curDir 

if [ -s $final ]; then
  echo job has already been done. if you want to do again, delete ${curDir}/${final} and submit job again.
  exit 0;
fi

echo STAR_start=`date` 

STAR $option --runThreadN $thread --genomeDir $genome_dir --readFilesIn $samples $uncompress --outFileNamePrefix ${sampleName}_ $output_format $rgline 

samtools index $final

samtools flagstat $final > ${final}.stat 

echo finished=`date` 

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
  my $sort_by_coordinate = get_option_value( $config->{$section}{sort_by_coordinate}, 0 );

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my @resultFiles = ();
    my $final = $sort_by_coordinate ? $sampleName . "_Aligned.sortedByCoord.out.bam" : $sampleName . "_Aligned.out.bam";
    push( @resultFiles, "${resultDir}/${sampleName}/${final}" );
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
