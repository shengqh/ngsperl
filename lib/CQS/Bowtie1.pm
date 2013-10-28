#!/usr/bin/perl
package CQS::Bowtie1;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::AbstractBowtie;

our @ISA = qw(CQS::AbstractBowtie);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "Bowtie1";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $bowtie1_index = $config->{$section}{bowtie1_index} or die "define ${section}::bowtie1_index first";
  my $samonly = $config->{$section}{samonly};
  if ( !defined $samonly ) {
    $samonly = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}_bt1.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . " \n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $samFile     = $sampleName . ".sam";
    my $bamFile     = $sampleName . ".bam";

    my $indent = "";
    my $tag    = "--sam-RG ID:$sampleName --sam-RG LB:$sampleName --sam-RG SM:$sampleName --sam-RG PL:ILLUMINA";

    my $fastqs = join( ',', @sampleFiles );
    my $bowtie1_aln_command = "bowtie $option -S $bowtie1_index $tag $fastqs -S $samFile";

    my $pbsName = "${sampleName}_bt1.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";
    my $curDir  = create_directory_or_die( $resultDir . "/$sampleName" );
    my $log     = "${logDir}/${sampleName}_bt1.log";

    print SH "\$MYCMD ./$pbsName \n";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

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
";
    }
    else {
      print OUT "
if [ -s $bamFile ]; then
  echo job has already been done. if you want to do again, delete $bamFile and submit job again.
  exit 0
fi

$bowtie1_aln_command

if [ -s $samFile ]; then
  samtools view -S -b $samFile | samtools sort - $sampleName
  samtools index $bamFile 
  samtools flagstat $bamFile > ${bamFile}.stat 
  rm $samFile
fi
";
    }

    print OUT "
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

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}
1;
