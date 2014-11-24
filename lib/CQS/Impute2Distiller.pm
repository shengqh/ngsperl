#!/usr/bin/perl
package CQS::Impute2Distiller;

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
  $self->{_name}   = "Impute2Distiller";
  $self->{_suffix} = "_id";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqsFile = get_cqstools( $config, $section, 1 );
  my $targetSnpFile = get_param_file( $config->{$section}{fasta_file}, "target_snp_file", 1 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @inputFiles   = @{ $rawFiles{$sampleName} };
    my $inputFileStr = join( ",", @inputFiles );

    my $genFile = $sampleName . ".gen";

    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir

if [ -s $genFile ]; then
  echo job has already been done. if you want to do again, delete $genFile and submit job again.
  exit 0
fi

echo Impute2Distiller=`date` 

mono-sgen $cqsFile impute2_distiller $option -i $inputFileStr -t $targetSnpFile -o $genFile

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

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $curDir = $resultDir . "/$sampleName";

    my @resultFiles = ();
    my $countFile   = "${curDir}/${sampleName}.gen";
    push( @resultFiles, $countFile );
    push( @resultFiles, "${curDir}/${sampleName}.info" );
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
