#!/usr/bin/perl
package CNV::CNVnator;

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
  $self->{_name}   = "CNV::CNVnator";
  $self->{_suffix} = "_cnvnator";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $binsize        = $config->{$section}{binsize}        or die "define ${section}::binsize first";
  my $chromosome_dir = $config->{$section}{chromosome_dir} or die "define ${section}::chromosome_dir first";

  my $genome    = $config->{$section}{genome};
  my $genomestr = "";

  if ( defined $genome ) {
    $genomestr = "-genome " . $genome;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $bamFile = $sampleFiles[0];

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    my $curDir   = create_directory_or_die( $resultDir . "/$sampleName" );
    my $rootFile = $sampleName . ".root";
    my $callFile = $sampleName . ".call";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $curDir

if [ -s $callFile ]; then
  echo job has already been done. if you want to do again, delete $callFile and submit job again.
else
  if [ ! -s $rootFile ]; then
    echo \"EXTRACTING READ MAPPING FROM BAM/SAM FILES =\" `date`
    cnvnator $genomestr -unique -root $rootFile -tree $bamFile 
  fi

  echo \"GENERATING HISTOGRAM =\" `date`
  cnvnator $genomestr -root $rootFile -d $chromosome_dir -his $binsize 

  echo \"CALCULATING STATISTICS =\" `date`
  cnvnator -root $rootFile -stat $binsize 

  echo \"RD SIGNAL PARTITIONING =\" `date`
  cnvnator -root $rootFile -partition $binsize 

  echo \"CNV CALLING =\" `date`
  cnvnator -root $rootFile -call $binsize > $callFile 
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

  print "!!!shell file $shfile created, you can run this shell file to submit all cnvnator tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( sort keys %rawFiles ) {
    my $curDir      = create_directory_or_die( $resultDir . "/$sampleName" );
    my $rootFile    = $sampleName . ".root";
    my $callFile    = $sampleName . ".call";
    my @resultFiles = ();
    push( @resultFiles, $curDir . "/" . $callFile );
    push( @resultFiles, $curDir . "/" . $rootFile );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }

  return $result;
}

1;
