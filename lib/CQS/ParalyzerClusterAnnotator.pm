#!/usr/bin/perl
package CQS::CQSMappedCount;

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
  $self->{_name} = "CQSMappedCount";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqsFile   = get_param_file( $config->{$section}{cqs_tools},  "cqs_tools",  1 );
  my $samtools  = get_param_file( $config->{$section}{samtools},   "samtools",   1 );
  my $gffFile   = get_param_file( $config->{$section}{gff_file},   "gff_file",   1 );
  my $fastaFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 0 );

  my $min_overlap = $config->{$section}{min_overlap};
  if ( defined $min_overlap ) {
    $option = $option . " --min_overlap $min_overlap";
  }

  if ( defined $fastaFile ) {
    $option = $option . " -f $fastaFile";
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my %seqCountFiles = ();
  if ( defined $config->{$section}{"seqcount"} || defined $config->{$section}{"seqcount_ref"} ) {
    %seqCountFiles = %{ get_raw_files( $config, $section, "seqcount" ) };
  }

  my %fastqFiles = ();
  if ( defined $config->{$section}{"fastq_files"} || defined $config->{$section}{"fastq_files_ref"} ) {
    %fastqFiles = %{ get_raw_files( $config, $section, "fastq_files" ) };
  }

  my $shfile = $pbsDir . "/${task_name}_ct.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @bamFiles  = @{ $rawFiles{$sampleName} };
    my $bamFile   = $bamFiles[0];
    my $fileName  = basename($bamFile);
    my $countFile = $fileName . ".count";

    my $seqcountFile = "";
    if ( defined $seqCountFiles{$sampleName} ) {
      my @seqcounts = @{ $seqCountFiles{$sampleName} };
      my $seqcount  = $seqcounts[0];
      $seqcountFile = " -c $seqcount";
    }

    my $fastqFile = "";
    if ( defined $fastqFiles{$sampleName} ) {
      my @files = @{ $fastqFiles{$sampleName} };
      my $file  = $files[0];
      $fastqFile = " -q $file";
    }

    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

    my $pbsName = "${sampleName}_ct.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${sampleName}_ct.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir

if [ -s $countFile ]; then
  echo job has already been done. if you want to do again, delete $countFile and submit job again.
  exit 0
fi

echo CQSMappedCount=`date` 

mono-sgen $cqsFile mapped_count $option --samtools $samtools -i $bamFile -g $gffFile -o $countFile $seqcountFile $fastqFile

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

  print "!!!shell file $shfile created, you can run this shell file to submit all trna_count tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $fasta_format = $config->{$section}{fasta_format};
  if ( !defined $fasta_format ) {
    $fasta_format = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $curDir = $resultDir . "/$sampleName";

    my @bamFiles = @{ $rawFiles{$sampleName} };
    my $bamFile  = $bamFiles[0];
    my $fileName = basename($bamFile);

    my @resultFiles = ();
    my $countFile   = "${curDir}/${fileName}.count";
    push( @resultFiles, $countFile );
    push( @resultFiles, "${countFile}.mapped.xml" );
    push( @resultFiles, "${curDir}/${fileName}.info" );

    my $unmapped;
    if ($fasta_format) {
      $unmapped = change_extension( $countFile, ".unmapped.fasta" );
    }
    else {
      $unmapped = change_extension( $countFile, ".unmapped.fastq" );
    }
    push( @resultFiles, $unmapped );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
