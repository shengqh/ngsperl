#!/usr/bin/perl
package Imputation::Impute2;

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
  $self->{_name}   = "Imputation::Impute2";
  $self->{_suffix} = "_imp";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  my %mapFiles = %{ get_raw_files( $config, $section, "map_file" ) };
  my %haploFiles = %{ get_raw_files( $config, $section, "haplo_file" ) };
  
  my $maxChromosomeLength = get_option($config, $section, "max_chromosome_length");
  my $interval = get_option($config, $section, "interval");

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $sample      = $sampleFiles[0];
    
    my @mapFiles = @{ $mapFiles{$sampleName} };
    my $map = $mapFiles[0];
    
    my @haploFiles = @{ $haploFiles{$sampleName} };
    my $haploFile = $haploFiles[0];
    my $legendFile = change_extension($haploFile, ".legend");

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

    my $gen_file  = "${sampleName}.gen";
    my $cat_command ="cat ";

    open( OUT, ">$pbsFile" ) or die $!;

    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir 

if [ -s $gen_file ]; then
  echo job has already been done. if you want to do again, delete $curDir/$gen_file and submit job again.
  exit 0;
fi

echo impute2_start=`date` 
";

    my $start = 1;
    while($start < $maxChromosomeLength){
      my $end = $start + $interval - 1;
      my $tmpFile = $sampleName . "_" . $start . ".tmp";
      $cat_command = $cat_command . $tmpFile . " ";
      print OUT "impute2 -known_haps_g $sample -m $map -int $start $end -h $haploFile -l $legendFile -o Pietenpol_p53.22.gen \n";
      $start = $end + 1;
    }

    print OUT "

$cat_command > $gen_file
rm *.tmp
echo finished=`date` 

";
    close(OUT);

    print SH "\$MYCMD ./$pbsName \n";
    print "$pbsFile created\n";
  }
  print SH "exit 0\n";
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my @resultFiles = ();

    push( @resultFiles, "${resultDir}/${sampleName}/${sampleName}.gen" );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
