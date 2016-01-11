#!/usr/bin/perl
package ATACseq::CleanBam;

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
  $self->{_name} = "ATACseq::CleanBam";
  $self->{_suffix} = "_cb";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $picard_jar = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1 );
  my $remove_chromosome = get_option($config, $section, "remove_chromosome"); 

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $sampleFile  = $sampleFiles[0];

    my $noChrFile = $sampleName . ".noChr" . $remove_chromosome. ".bam";
    my $finalFile = $sampleName . ".noChr" . $remove_chromosome. ".rmdup.bam";

    my $pbsFile = $self->pbsfile($pbsDir, $sampleName);
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

if [ -s $finalFile ]; then
  echo job has already been done. if you want to do again, delete $finalFile and submit job again.
  exit 0
fi

if [ ! -s $noChrFile ]; then
  echo RemoveChromosome=`date` 
  samtools idxstats $sampleFile | cut -f 1 | grep -v $remove_chromosome | xargs samtools view -b $sampleFile > $noChrFile 
fi

if [[ -s $noChrFile && ! -s $finalFile ]]; then
  echo RemoveDuplicate=`date` 
  java $option -jar $picard_jar MarkDuplicates I=$noChrFile O=$finalFile ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT M=${finalFile}.metrics
fi

if [[ -s $finalFile && ! -s ${finalFile}.bai ]]; then
  echo BamIndex=`date` 
  samtools index $finalFile
  samtools flagstat $finalFile > ${finalFile}.stat
  rm $noChrFile
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

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  my $remove_chromosome = get_option($config, $section, "remove_chromosome"); 

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $sortedFile   = $sampleName . ".noChr" . $remove_chromosome. ".rmdup.bam";

    my @resultFiles = ();
    push( @resultFiles, "${resultDir}/${sampleName}/${sortedFile}" );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
