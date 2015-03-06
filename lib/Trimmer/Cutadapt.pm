#!/usr/bin/perl
package Trimmer::Cutadapt;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "Cutadapt";
  $self->{_suffix} = "_cut";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $adapter   = get_option( $config, $section, "adapter" );
  my $extension = get_option( $config, $section, "extension" );
  my $gzipped   = get_option( $config, $section, "gzipped", 1 );

  if ( $gzipped && $extension =~ /\.gz$/ ) {

    #print $extension . "\n";
    $extension =~ s/\.gz$//g;

    #print $extension . "\n";
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  my $shortLimited = $option =~ /-m\s+\d+/;
  my $longLimited  = $option =~ /-M\s+\d+/;

  for my $sampleName ( sort keys %rawFiles ) {

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir
";

    my @sampleFiles = @{ $rawFiles{$sampleName} };
    if ( scalar(@sampleFiles) == 1 ) {
      my $finalName      = $sampleName . $extension;
      my $finalShortName = $finalName . ".short";
      my $finalLongName  = $finalName . ".long";

      my $finalFile      = $gzipped ? "${finalName}.gz"      : $finalName;
      my $finalShortFile = $gzipped ? "${finalShortName}.gz" : $finalShortName;
      my $finalLongFile  = $gzipped ? "${finalLongName}.gz"  : $finalLongName;

      print OUT "
if [ -s $finalFile ];then
  echo job has already been done. if you want to do again, delete ${resultDir}/$finalFile and submit job again.
  exit 0;
fi

";
      print OUT "cutadapt $option -a $adapter -o $finalFile ";
      if ($shortLimited) {
        print OUT " --too-short-output=$finalShortFile";
      }
      if ($longLimited) {
        print OUT " --too-long-output=$finalLongFile";
      }
      print OUT " $sampleFiles[0] \n";
    }
    else {

      #pair-end data
      my $read1file = $sampleFiles[0];
      my $read2file = $sampleFiles[1];
      my $read1name = $sampleName . ".1.fastq.gz";
      my $read2name = $sampleName . ".2.fastq.gz";

      print OUT "
if [ -s $read1name ];then
  echo job has already been done. if you want to do again, delete ${resultDir}/$read1name and submit job again.
  exit 0;
fi

";
      if ( $shortLimited || $longLimited ) {
        my $temp1name = $sampleName . ".1.tmp.fastq";
        my $temp2name = $sampleName . ".2.tmp.fastq";

        #https://cutadapt.readthedocs.org/en/stable/guide.html#illumina-truseq
        print OUT "cutadapt $option -a $adapter -o $temp1name -p $temp2name $read1file $read2file \n";
        print OUT "cutadapt $option -a $adapter -o $read2name -p $read1name $temp2name $temp1name \n";
        print OUT "rm $temp2name $temp1name \n";
      }
      else {
        print OUT "cutadapt $option -a $adapter -o $read1name $read1file \n";
        print OUT "cutadapt $option -a $adapter -o $read2name $read2file \n";
      }
    }
    print OUT "
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

  my $extension = $config->{$section}{extension} or die "define ${section}::extension first";
  my $gzipped = get_option( $config, $section, "gzipped", 1 );

  if ( $gzipped && $extension =~ /\.gz$/ ) {
    $extension =~ s/\.gz$//g;
  }

  my $shortLimited = $option =~ /-m\s+\d+/;
  my $longLimited  = $option =~ /-M\s+\d+/;

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result      = {};

  for my $sampleName ( keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my @resultFiles = ();
    if ( scalar(@sampleFiles) == 1 ) {
      my $finalName      = $sampleName . $extension;
      my $finalShortName = $finalName . ".short";
      my $finalLongName  = $finalName . ".long";

      my $finalFile      = $gzipped ? "${finalName}.gz"      : $finalName;
      my $finalShortFile = $gzipped ? "${finalShortName}.gz" : $finalShortName;
      my $finalLongFile  = $gzipped ? "${finalLongName}.gz"  : $finalLongName;

      push( @resultFiles, $resultDir . "/" . $finalFile );

      if ($shortLimited) {
        push( @resultFiles, $resultDir . "/" . $finalShortFile );
      }
      if ($longLimited) {
        push( @resultFiles, $resultDir . "/" . $finalLongFile );
      }
    }
    else {

      #pair-end data
      my $read1name = $sampleName . ".1.fastq.gz";
      my $read2name = $sampleName . ".2.fastq.gz";

      push( @resultFiles, $resultDir . "/" . $read1name );
      push( @resultFiles, $resultDir . "/" . $read2name );
    }
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
