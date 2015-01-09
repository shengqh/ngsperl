#!/usr/bin/perl
package Picard::MarkDuplicates;

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
  $self->{_name} = "Picard::MarkDuplicates";
  $self->{_suffix} = "_rd";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $markDuplicates_jar = get_param_file( $config->{$section}{markDuplicates_jar}, "markDuplicates_jar", 1 );
  my $thread_count = $config->{$section}{thread_count};

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $sampleFile  = $sampleFiles[0];

    my $rmdupFile    = $sampleName . ".rmdup.bam";
    my $sortedPrefix = $sampleName . ".rmdup_sorted";
    my $sortedFile   = $sortedPrefix . ".bam";

    my $pbsFile = $self->pbsfile($pbsDir, $sampleName);
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    my $curDir  = create_directory_or_die( $resultDir . "/$sampleName" );

    print SH "\$MYCMD ./$pbsName \n";

    my $cluster = get_cluster( $config, $section );
    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $curDir

if [ -s $sortedFile ]; then
  echo job has already been done. if you want to do again, delete $sortedFile and submit job again.
  exit 0
fi

if [ ! -s $rmdupFile ]; then
  echo RemoveDuplicate=`date` 
  java $option -jar $markDuplicates_jar I=$sampleFile O=$rmdupFile M=${rmdupFile}.matrix VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true
fi

if [[ -s $rmdupFile && ! -s $sortedFile ]]; then
  echo BamSort=`date` 
  samtools sort $rmdupFile $sortedPrefix 
fi

if [[ -s $sortedFile && ! -s ${sortedFile}.bai ]]; then
  echo BamIndex=`date` 
  samtools index $sortedFile
  samtools flagstat $sortedFile > ${sortedFile}.stat
  rm $rmdupFile
fi
  
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

  print "!!!shell file $shfile created, you can run this shell file to submit all Picard::MarkDuplicates tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $sortedFile   = $sampleName . ".rmdup_sorted.bam";

    my @resultFiles = ();
    push( @resultFiles, "${resultDir}/${sampleName}/${sortedFile}" );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
