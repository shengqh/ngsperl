#!/usr/bin/perl
package CQS::SmallRNACount;

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
  $self->{_name}   = "SmallRNACount";
  $self->{_suffix} = "_sc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $cqsFile = get_cqstools( $config, $section, 1 );
  my $coordinate_file = get_param_file( $config->{$section}{coordinate_file}, "coordinate_file", 1 );
  my $fastaFile       = get_param_file( $config->{$section}{fasta_file},      "fasta_file",      0 );

  if ( defined $fastaFile ) {
    $option = $option . " -f $fastaFile";
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my %seqCountFiles = ();
  if ( has_raw_files($config, $section, "seqcount")) {
    %seqCountFiles = %{ get_raw_files( $config, $section, "seqcount" ) };
  }

  my %fastqFiles = ();
  if ( has_raw_files($config, $section, "fastq_files")) {
    %fastqFiles = %{ get_raw_files( $config, $section, "fastq_files" ) };
  }

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @bamFiles  = @{ $rawFiles{$sampleName} };
    my $bamFile   = $bamFiles[0];
    my $countFile = $sampleName . ".count";

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

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $curDir

if [ -s $countFile ]; then
  echo job has already been done. if you want to do again, delete $countFile and submit job again.
  exit 0
fi

echo SmallRNACount=`date` 

mono $cqsFile smallrna_count $option -i $bamFile -g $coordinate_file $seqcountFile $fastqFile -o $countFile

echo finished=`date`

exit 0 
";

    close OUT;

    print "$pbsFile created \n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  my $unmapped_fastq = $option =~ /--unmapped_fastq/;

  my %seqCountFiles = ();
  if ( defined $config->{$section}{"seqcount"} || defined $config->{$section}{"seqcount_ref"} ) {
    %seqCountFiles = %{ get_raw_files( $config, $section, "seqcount" ) };
  }

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $curDir = $resultDir . "/$sampleName";

    my @resultFiles = ();
    my $countFile   = "${curDir}/${sampleName}.count";
    my $tRNAPositionFile   = "${curDir}/${sampleName}.tRNA.position";
    push( @resultFiles, $countFile );
    push( @resultFiles, $tRNAPositionFile );
    push( @resultFiles, "${countFile}.mapped.xml" );
    push( @resultFiles, "${curDir}/${sampleName}.info" );

    if ($unmapped_fastq) {
      my $unmapped = change_extension( $countFile, ".unmapped.fastq.gz" );
      push( @resultFiles, $unmapped );
      if ( defined $seqCountFiles{$sampleName} ) {
        push( @resultFiles, $unmapped . ".dupcount" );
      }
    }

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
