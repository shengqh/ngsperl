#!/usr/bin/perl
package CQS::CNV;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::DNASeq;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::ClassFactory;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(cnvnator conifer cnmops freec)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub cnvnator {
  my ( $config, $section ) = @_;
  my $obj = instantiate("CNV::CNVnator");
  $obj->perform( $config, $section );
}

sub conifer {
  my ( $config, $section ) = @_;
  my $obj = instantiate("CNV::Conifer");
  $obj->perform( $config, $section );
}

sub cnmops {
  my ( $config, $section ) = @_;
  my $obj = instantiate("CNV::cnMops");
  $obj->perform( $config, $section );
}

sub freec {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );
  my $chrLenFile             = $config->{$section}{chrLenFile}             or die "define ${section}::chrLenFile first";
  my $ploidy                 = $config->{$section}{ploidy}                 or die "define ${section}::ploidy first";
  my $coefficientOfVariation = $config->{$section}{coefficientOfVariation} or die "define ${section}::coefficientOfVariation first";
  my $chrFiles               = $config->{$section}{chrFiles}               or die "define ${section}::chrFiles first";
  my $inputFormat            = $config->{$section}{inputFormat}            or die "define ${section}::inputFormat first";
  my $mateOrientation        = $config->{$section}{mateOrientation}        or die "define ${section}::mateOrientation first";
  my $bedfile = $config->{$section}{bedfile};

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $bamFile = $sampleFiles[0];

    my $pbsName = "freec_${sampleName}.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/freec_${sampleName}.log";
    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );
    my $configName = "${sampleName}.conf";
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

echo CNV_CALLING=`date`
freec -conf $configName 

echo finished=`date`
";
    close OUT;

    print "$pbsFile created\n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all freec tasks.\n";

  #`qsub $pbsFile`;
}

1;
