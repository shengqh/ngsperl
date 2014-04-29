#!/usr/bin/perl
package CNV::Freec;

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
  $self->{_name}   = "CNV::Freec";
  $self->{_suffix} = "_fc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $chrLenFile             = $config->{$section}{chrLenFile}             or die "define ${section}::chrLenFile first";
  my $chrFiles               = $config->{$section}{chrFiles}               or die "define ${section}::chrFiles first, it is the path of directory contains chromosome fasta files";
  my $ploidy                 = $config->{$section}{ploidy}                 or die "define ${section}::ploidy first, such as 2";
  my $coefficientOfVariation = $config->{$section}{coefficientOfVariation} or die "define ${section}::coefficientOfVariation first, such as 0.05";
  my $inputFormat            = $config->{$section}{inputFormat}            or die "define ${section}::inputFormat first";
  my $mateOrientation        = $config->{$section}{mateOrientation}
    or die "define ${section}::mateOrientation first";    #0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)
  my $maxThreads = $config->{$section}{maxThreads} or die "define ${section}::maxThreads first";

  my $bedfile = $config->{$section}{bedfile};

  my $rawFiles = get_raw_files( $config, $section );
  print(%{$rawFiles}, "\n");
  my $groups = get_raw_files( $config, $section, "groups" );

  my $shfile = $pbsDir . "/${task_name}.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $groupName ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$groupName} };
    print(@samples, "\n");

    my $bamFile;
    my $bamFile2;
    if ( scalar(@samples) > 1 ) {
      $bamFile2 = $rawFiles->{ $samples[1] }[1];    #control
      $bamFile  = $rawFiles->{ $samples[2] }[1];    #sample
    }
    else {
      my $bamFile = $rawFiles->{ $samples[1] }[1];    #sample
    }

    my $pbsFile = $self->pbsfile( $pbsDir, $groupName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $groupName );

    print SH "\$MYCMD ./$pbsName \n";
    my $curDir     = create_directory_or_die( $resultDir . "/$groupName" );
    my $configName = "${groupName}.conf";
    my $configFile = ${curDir} . "/$configName";

    open( CON, ">$configFile" ) or die $!;

    print CON "[general] 
chrLenFile = $chrLenFile 
ploidy = $ploidy 
coefficientOfVariation = $coefficientOfVariation 
chrFiles = $chrFiles 

[sample]
mateFile = $bamFile 
inputFormat  = $inputFormat 
mateOrientation = $mateOrientation 

";
    if ( defined $bamFile2 ) {
      print CON "[control] 
mateFile = $bamFile2 
inputFormat  = $inputFormat 
mateOrientation = $mateOrientation 

";

    }

    if ( defined $bedfile ) {
      print CON "[target] \n\n";
      print CON "captureRegions = $bedfile \n\n";
    }

    close(CON);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir

echo freec_start=`date`
freec -conf $configName 

echo freec_end=`date`
";
    close OUT;

    print "$pbsFile created\n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all freec tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = { $task_name => [ $resultDir . "/${task_name}.call" ] };

  return $result;
}

1;
