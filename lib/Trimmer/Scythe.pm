#!/usr/bin/perl
package Trimmer::Scythe;

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
  $self->{_name}   = "Trimmer::Scythe";
  $self->{_suffix} = "_scy";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $faFile = get_param_file( $config->{$section}{adapter_file}, "adapter_file", 1 );
  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

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

";
    for my $sampleFile (@sampleFiles) {
      my $sname = basename($sampleFile);
      if($sname =~ ".gz\$"){
        $sname =  substr $sname, 1, -3;
      }
      my $trimFile = change_extension( $sname, "_trim.fastq" );
      my $finalFile = $trimFile . ".gz";

      print OUT "if [ ! -s $finalFile ]; then
  scythe $option -a $faFile $sampleFile > $trimFile
  gzip $trimFile
fi
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

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name}  . " tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my @resultFiles = ();

    for my $sampleFile (@sampleFiles) {
      my $sname = basename($sampleFile);
      my $trimFile = change_extension( $sname, "_trim.fastq" );
      my $finalFile = $trimFile . ".gz";
      push( @resultFiles, "${resultDir}/${finalFile}" );
    }
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
