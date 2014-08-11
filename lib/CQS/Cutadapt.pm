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
  $self->{_name}   = "Cutadapt";
  $self->{_suffix} = "_cut";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $adapt     = get_option( $config, $section, "adaptor" );
  my $extension = get_option( $config, $section, "extension" );
  my $gzipped   = get_option( $config, $section, "gzipped", 1 );
  
  if($gzipped && $extension =~ /\.gz$/){
    print $extension . "\n"; 
    $extension =~ s/\.gz$//g;
    print $extension . "\n"; 
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  my $shortLimited = $option =~ /-m\s+\d+/;
  my $longLimited  = $option =~ /-M\s+\d+/;

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $finalName      = $sampleName . $extension;
    my $finalShortName = $finalName . ".short";
    my $finalLongName  = $finalName . ".long";

    my $finalFile      = $gzipped ? "${finalName}.gz"      : $finalName;
    my $finalShortFile = $gzipped ? "${finalShortName}.gz" : $finalShortName;
    my $finalLongFile  = $gzipped ? "${finalLongName}.gz"  : $finalLongName;

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir

if [ -s $finalFile ];then
  echo job has already been done. if you want to do again, delete ${resultDir}/$finalFile and submit job again.
  exit 0;
fi

";

    if ( scalar(@sampleFiles) == 1 ) {
      print OUT "cutadapt $sampleFiles[0] $option -a $adapt -o $finalName";
      if ($shortLimited) {
        print OUT " --too-short-output=$finalShortName";
      }
      if ($longLimited) {
        print OUT " --too-long-output=$finalLongName";
      }
      print OUT "\n";
    }
    else {
      my $outputFiles = "";
      my $shortFiles  = "";
      my $longFiles   = "";
      for my $sampleFile (@sampleFiles) {
        my $fileName   = basename($sampleFile);
        my $outputFile = change_extension( $fileName, $extension );
        my $shortFile  = $outputFile . ".short";
        my $longFile   = $outputFile . ".long";
        $outputFiles = $outputFiles . " " . $outputFile;
        $shortFiles  = $shortFiles . " " . $shortFile;
        $longFiles   = $longFiles . " " . $longFile;
        print OUT "cutadapt $sampleFile $option -a $adapt -o $outputFile";

        if ($shortLimited) {
          print OUT " --too-short-output=$shortFile";
        }
        if ($longLimited) {
          print OUT " --too-long-output=$longFile";
        }
        print OUT "\n";
        print OUT "
cat $outputFiles > $finalName
rm $outputFiles
";

        if ($shortLimited) {
          print OUT "
cat $shortFiles > $finalShortName
rm $shortFiles
";
        }
        if ($longLimited) {
          print OUT "
cat $longFiles > $finalLongName
rm $longFiles
";
        }
      }
    }
    if ($gzipped) {
      print OUT "
gzip $finalName
";
      if ($shortLimited) {
        print OUT "
gzip $finalShortName
";
      }
      if ($longLimited) {
        print OUT "
gzip $finalLongName
";
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

  print "!!!shell file $shfile created, you can run this shell file to submit all bwa tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $extension = $config->{$section}{extension} or die "define ${section}::extension first";
  my $gzipped = get_option( $config, $section, "gzipped", 1 );
  my $shortLimited = $option =~ /-m\s+\d+/;
  my $longLimited  = $option =~ /-M\s+\d+/;

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $finalName = $sampleName . $extension;

    my $finalFile = $gzipped ? "${finalName}.gz" : $finalName;

    my @resultFiles = ();
    push( @resultFiles, $resultDir . "/" . $finalFile );

    if ($shortLimited) {
      my $finalShortName = $finalName . ".short";
      my $finalShortFile = $gzipped ? "${finalShortName}.gz" : $finalShortName;
      push( @resultFiles, $resultDir . "/" . $finalShortFile );
    }

    if ($longLimited) {
      my $finalLongName = $finalName . ".long";
      my $finalLongFile = $gzipped ? "${finalLongName}.gz" : $finalLongName;
      push( @resultFiles, $resultDir . "/" . $finalLongFile );
    }

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;