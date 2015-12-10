#!/usr/bin/perl
package Alignment::Bowtie1;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use Alignment::AbstractBowtie;

our @ISA = qw(Alignment::AbstractBowtie);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "Bowtie1";
  $self->{_suffix} = "_bt1";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $bowtie1_index = $config->{$section}{bowtie1_index} or die "define ${section}::bowtie1_index first";
  my $samformat  = get_option( $config, $section, "samformat",  1 );
  my $mappedonly = get_option( $config, $section, "mappedonly", 0 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . " \n";

  if ($samformat) {
    my $samonly = get_option( $config, $section, "samonly", 0 );
    my $sortbam = get_option( $config, $section, "sortbam", 1 );

    for my $sampleName ( sort keys %rawFiles ) {
      my @sampleFiles = @{ $rawFiles{$sampleName} };
      my $samFile     = $sampleName . ".sam";

      my $bowtiesam        = $samFile;
      my $mappedonlycmd    = "";
      my $mappedonlyoption = "";
      if ($mappedonly) {
        $bowtiesam     = $sampleName . ".all.sam";
        $mappedonlycmd = "
if [ -s $bowtiesam ]; then
  samtools view -F 4 $bowtiesam > $samFile
  rm $bowtiesam
fi";
        $mappedonlyoption = "-F 4";
      }

      my $bamFile = $sampleName . ".bam";

      my $indent = "";
      my $tag    = "--sam-RG ID:$sampleName --sam-RG LB:$sampleName --sam-RG SM:$sampleName --sam-RG PL:ILLUMINA";

      my $bowtie1_aln_command;
      if ( $sampleFiles[0] =~ /.gz$/ ) {
        if ( scalar(@sampleFiles) == 1 ) {
          $bowtie1_aln_command = "zcat $sampleFiles[0] | bowtie $option -S $tag $bowtie1_index - $bowtiesam";
        }
        else {
          my $f1 = $sampleFiles[0];
          my $f2 = $sampleFiles[1];

          $bowtie1_aln_command = "
mkfifo ${f1}.fifo
zcat $f1 > ${f1}.fifo &

mkfifo ${f2}.fifo
zcat $f2 > ${f2}.fifo &
        
bowtie $option -S $tag $bowtie1_index ${f1}.fifo,${f2}.fifo $bowtiesam
 
rm ${f1}.fifo
rm ${f2}.fifo";
        }
      }
      else {
        my $fastqs = join( ',', @sampleFiles );
        $bowtie1_aln_command = "bowtie $option -S $tag $bowtie1_index $fastqs $bowtiesam";
      }

      my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
      my $pbsName = basename($pbsFile);
      my $log     = $self->logfile( $logDir, $sampleName );
      my $curDir  = create_directory_or_die( $resultDir . "/$sampleName" );

      print SH "\$MYCMD ./$pbsName \n";

      my $log_desc = $cluster->get_log_desc($log);

      open( OUT, ">$pbsFile" ) or die $!;
      print OUT "$pbsDesc
$log_desc

$path_file

cd $curDir
";
      if ($samonly) {
        print OUT "
if [ -s $samFile ]; then
  echo job has already been done. if you want to do again, delete $samFile and submit job again.
  exit 0
fi

$bowtie1_aln_command

$mappedonlycmd
";
      }
      else {
        print OUT "
if [ -s $bamFile ]; then
  echo job has already been done. if you want to do again, delete $bamFile and submit job again.
  exit 0
fi

$bowtie1_aln_command

if [ -s $bowtiesam ]; then
";
        if ($sortbam) {
          print OUT "  samtools view -Shu $mappedonlyoption $bowtiesam | samtools sort -o $bamFile -
  if [ -s $bamFile ]; then
    samtools index $bamFile 
    samtools flagstat $bamFile > ${bamFile}.stat
";
        }
        else {
          print OUT "samtools view -S $mappedonlyoption -b $bowtiesam > ${sampleName}.bam
  if [ -s $bamFile ]; then
";
        }
        print OUT "    rm $bowtiesam
  fi
fi
";
      }

      print OUT "
echo finished=`date`

exit 0;
";

      close OUT;

      print "$pbsFile created\n";
    }
    close(SH);

    if ( is_linux() ) {
      chmod 0755, $shfile;
    }
  }
  else {
    for my $sampleName ( sort keys %rawFiles ) {
      my @sampleFiles = @{ $rawFiles{$sampleName} };
      my $finalFile   = $sampleName . ".out";

      my $indent = "";

      my $bowtie1_aln_command;
      if ( $sampleFiles[0] =~ /.gz$/ ) {
        if ( scalar(@sampleFiles) == 1 ) {
          $bowtie1_aln_command = "zcat $sampleFiles[0] | bowtie $option $bowtie1_index - $finalFile";
        }
        else {
          my $f1 = $sampleFiles[0];
          my $f2 = $sampleFiles[1];

          $bowtie1_aln_command = "
mkfifo ${f1}.fifo
zcat $f1 > ${f1}.fifo &

mkfifo ${f2}.fifo
zcat $f2 > ${f2}.fifo &
        
bowtie $option $bowtie1_index ${f1}.fifo,${f2}.fifo $finalFile
 
rm ${f1}.fifo
rm ${f2}.fifo
";
        }
      }
      else {
        my $fastqs = join( ',', @sampleFiles );
        $bowtie1_aln_command = "bowtie $option $bowtie1_index $fastqs $finalFile ";
      }

      my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
      my $pbsName = basename($pbsFile);
      my $log     = $self->logfile( $logDir, $sampleName );

      my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

      print SH "\$MYCMD ./$pbsName \n";

      open( OUT, ">$pbsFile" ) or die $!;
      print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir

if [ -s $finalFile ]; then
  echo job has already been done. if you want to do again, delete $finalFile and submit job again.
  exit 0
fi

$bowtie1_aln_command

echo finished=`date`

exit 0;
";

      close OUT;

      print "$pbsFile created\n";
    }
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

1;
