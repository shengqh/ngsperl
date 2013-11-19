#!/usr/bin/perl
package CQS::BWA;

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
  $self->{_name} = "BWA";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}_bwa.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $samFile     = $sampleName . ".sam";
    my $bamFile     = $sampleName . ".bam";
    my $tag         = get_bam_tag($sampleName);

    my $sampleFile1 = $sampleFiles[0];

    my $bwa_aln_command;
    if ( scalar(@sampleFiles) == 2 ) {
      my $sampleFile2 = $sampleFiles[1];
      $bwa_aln_command = "
if [ ! -s $samFile ]; then
  echo bwa_mem=`date` 
  bwa mem $option -r $tag $faFile $sampleFile1 $sampleFile2 > $samFile
fi";
    }
    else {
      $bwa_aln_command = "
if [ ! -s $samFile ]; then
  echo bwa_mem=`date` 
  bwa mem $option -r $tag $faFile $sampleFile1 > $samFile
fi";
    }

    my $pbsName = "${sampleName}_bwa.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";
    my $curDir  = create_directory_or_die( $resultDir . "/$sampleName" );
    my $log     = "${logDir}/${sampleName}_bwa.log";

    print SH "\$MYCMD ./$pbsName \n";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir

if [ -s $bamFile ]; then
  echo job has already been done. if you want to do again, delete $bamFile and submit job again.
  exit 0
fi

$bwa_aln_command

if [ -s $samFile ]; then
  samtools view -S -b $samFile | samtools sort - $sampleName
  samtools index $bamFile 
  samtools flagstat $bamFile > ${bamFile}.stat 
  rm $samFile
fi
  
echo finished=`date`

exit 1;
";

    close OUT;

    print "$pbsFile created\n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all BWA tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $bamFile     = "${resultDir}/${sampleName}/${sampleName}.bam";
    my @resultFiles = ();
    push( @resultFiles, $bamFile );
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
