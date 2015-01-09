#!/usr/bin/perl
package CNV::Conifer;

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
  $self->{_name}   = "CNV::Conifer";
  $self->{_suffix} = "_conifer";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $conifer = $config->{$section}{conifer} or die "define conifer program location first.\nconifer => \"location\"";
  my $bedfile = $config->{$section}{bedfile} or die "define $section:bedfile first";

  my $isbamsorted = $config->{$section}{isbamsorted};
  if ( !defined($isbamsorted) ) {
    $isbamsorted = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  create_directory_or_die( $resultDir . "/rpkm" );
  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    my $bamFile = $sampleFiles[0];

    if ( !$isbamsorted ) {
      ( $bamFile, my $bamSorted ) = get_sorted_bam($bamFile);

      #print $bamFile . "\n";
    }
    my $rpkm = "rpkm/" . $sampleName . ".rpkm";

    my $cluster = get_cluster( $config, $section );
    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

if [ ! -s $rpkm ]; then
  echo conifer_RPKM=`date`
  python $conifer rpkm --probes $bedfile --input $bamFile --output $rpkm 
fi

echo finished=`date`
";
    close OUT;

    print "$pbsFile created\n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  my $pbsFile   = $self->pbsfile( $pbsDir, $task_name . "_after_rpkm" );
  my $pbsName   = basename($pbsFile);
  my $hdf5File  = "${task_name}.hdf5";
  my $svalsFile = "${task_name}.svals";
  my $plotFile  = "${task_name}.png";
  my $sdFile    = "${task_name}.sd";
  my $callFile  = "${task_name}.call";

  my $log = "${logDir}/${task_name}_after_rpkm.log";
  create_directory_or_die( $resultDir . "/call_images" );

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir

#2 analysis
echo analyze=`date`
python $conifer analyze --probes $bedfile --rpkm_dir rpkm/ --output $hdf5File --svd 6 --write_svals $svalsFile --plot_scree $plotFile --write_sd $sdFile 

#3 call
echo call=`date`
python $conifer call --input $hdf5File --output $callFile 

#4 plot
python $conifer plotcalls --input $hdf5File --calls $callFile --output call_images 
";
  close OUT;

  print "$pbsFile created\n";

  print "!!!shell file $shfile created, you can run this shell file to submit all conifer rpkm tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result = { $task_name => [ $resultDir . "/${task_name}.call" ] };

  return $result;
}

sub pbsfiles {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %fqFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( sort keys %fqFiles ) {
    $result->{$sampleName} = $self->pbsfile( $pbsDir, $sampleName );
  }

  $result->{$task_name} = $self->pbsfile( $pbsDir, $task_name . "_after_rpkm" );

  return $result;
}

1;
