#!/usr/bin/perl
package Alignment::Bowtie2;

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
  $self->{_name}   = "Bowtie2";
  $self->{_suffix} = "_bt2";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $bowtie2_index = $config->{$section}{bowtie2_index} or die "define ${section}::bowtie2_index first";
  my $samonly = get_option($config, $section, "samonly", 0);

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct), "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $samFile     = $sampleName . ".sam";
    my $bamFile     = $sampleName . ".bam";

    my $indent = "";
    my $tag    = "--sam-RG ID:$sampleName --sam-RG LB:$sampleName --sam-RG SM:$sampleName --sam-RG PL:ILLUMINA";

    my $fastqs = join( ',', @sampleFiles );
    my $bowtie2_aln_command = "bowtie2 $option -x $bowtie2_index -U $fastqs $tag -S $samFile";

    my $index_command = get_index_command( $bamFile, $indent );
    my $stat_command = get_stat_command( $bamFile, $indent );

    my $pbsName = $self->pbsname($sampleName);
    my $pbsFile = $pbsDir . "/$pbsName";
    my $log     = $self->logfile( $logDir, $sampleName );

    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

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

$bowtie2_aln_command
";
    }
    else {
      print OUT "
if [ -s $bamFile ]; then
  echo job has already been done. if you want to do again, delete $bamFile and submit job again.
  exit 0
fi

$bowtie2_aln_command

if [ -s $samFile ]; then
  samtools view -S -b -F 256 $samFile | samtools sort - $sampleName
  if [ -s $bamFile ]; then
    samtools index $bamFile 
    samtools flagstat $bamFile > ${bamFile}.stat 
    rm $samFile
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

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

1;
