#!/usr/bin/perl
package CQS::Cutadapt;

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
  $self->{_name} = "Cutadapt";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $adapt     = $config->{$section}{adaptor}   or die "define ${section}::adaptor first";
  my $extension = $config->{$section}{extension} or die "define ${section}::extension first";
  
  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}_cutadapt.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $finalFile   = $sampleName . $extension;

    my $pbsName = "${sampleName}_cut.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${sampleName}_cut.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir

if [ -s $finalFile ];then
  echo job has already been done. if you want to do again, delete ${resultDir}/${finalFile} and submit job again.
  exit 1;
fi

";
    if ( scalar(@sampleFiles) == 1 ) {
      print OUT "cutadapt $sampleFiles[0] $option -a $adapt -o $finalFile \n";
    }
    else {
      my $outputFiles = "";
      for my $sampleFile (@sampleFiles) {
        my $fileName = basename($sampleFile);
        my $outputFile = change_extension( $fileName, $extension );
        $outputFiles = $outputFiles . " " . $outputFile;
        print OUT "cutadapt $sampleFile $option -a $adapt -o $outputFile \n";
      }

      print OUT "
cat $outputFiles > $finalFile
 
rm $outputFiles
";
    }
    print OUT "
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

  print "!!!shell file $shfile created, you can run this shell file to submit all bwa tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $extension = $config->{$section}{extension} or die "define ${section}::extension first";

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $finalFile   = $sampleName . $extension;
    my @resultFiles = ();
    push( @resultFiles, $resultDir . "/" . $finalFile );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
