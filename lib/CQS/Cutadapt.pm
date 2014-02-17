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
  $self->{_suffix} = "_cut";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $adapt     = get_option( $config, $section, "adaptor" );
  my $extension = get_option( $config, $section, "extension" );
  my $gzipped   = get_option( $config, $section, "gzipped", 0 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $finalName       = $sampleName . $extension;
    my $finalShortName  = $finalName . ".short";
    my $finalUntrimName = $finalName . ".untrimmed";

    my $finalFile       = $gzipped ? "${finalName}.gz"       : $finalName;
    my $finalShortFile  = $gzipped ? "${finalShortName}.gz"  : $finalShortName;
    my $finalUntrimFile = $gzipped ? "${finalUntrimName}.gz" : $finalUntrimName;

    my $pbsName = $self->pbsname($sampleName);
    my $pbsFile = $pbsDir . "/$pbsName";
    my $log     = $self->logname( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir

if [ -s $finalFile ];then
  echo job has already been done. if you want to do again, delete ${resultDir}/$finalFile and submit job again.
  exit 1;
fi

";
    if ( scalar(@sampleFiles) == 1 ) {
      print OUT "cutadapt $sampleFiles[0] $option -a $adapt -o $finalName --too-short-output=$finalShortName --untrimmed-output=$finalUntrimName \n";
    }
    else {
      my $outputFiles = "";
      my $shortFiles  = "";
      my $untrimFiles = "";
      for my $sampleFile (@sampleFiles) {
        my $fileName   = basename($sampleFile);
        my $outputFile = change_extension( $fileName, $extension );
        my $shortFile  = $outputFile . ".short";
        my $untrimFile = $outputFile . ".untrimed";
        $outputFiles = $outputFiles . " " . $outputFile;
        $shortFiles  = $shortFiles . " " . $shortFile;
        $untrimFiles = $untrimFiles . " " . $untrimFile;
        print OUT "cutadapt $sampleFile $option -a $adapt -o $outputFile --too-short-output=$shortFile --untrimmed-output=$untrimFile \n";
      }

      print OUT "
cat $outputFiles > $finalName
rm $outputFiles
cat $shortFiles > $finalShortName
rm $shortFiles
cat $untrimFiles > $finalUntrimName
rm $untrimFiles
";

    }
    if ($gzipped) {
      print OUT "
gzip $finalName
gzip $finalShortName
gzip $finalUntrimName
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
  my $gzipped = get_option( $config, $section, "gzipped", 1 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $finalName = $sampleName . $extension;

    my $finalFile = $gzipped ? "${finalName}.gz" : $finalName;

    my @resultFiles = ();
    push( @resultFiles, $resultDir . "/" . $finalFile );
    if ( $option !~ "-M" ) {
      my $finalShortName  = $finalName . ".short";
      my $finalUntrimName = $finalName . ".untrimmed";
      my $finalShortFile  = $gzipped ? "${finalShortName}.gz" : $finalShortName;
      my $finalUntrimFile = $gzipped ? "${finalUntrimName}.gz" : $finalUntrimName;
      push( @resultFiles, $resultDir . "/" . $finalShortFile );
      push( @resultFiles, $resultDir . "/" . $finalUntrimFile );
    }

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
