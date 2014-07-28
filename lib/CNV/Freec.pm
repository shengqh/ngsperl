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

  my $snpFile = parse_param_file( $config, $section, "SNPfile", 0 );
  if ( defined $snpFile ) {
    $self->performPileup( $config, $section );
  }
  else {
    $self->performBAM( $config, $section );
  }
}

sub performBAM {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $chrLenFile             = $config->{$section}{chrLenFile}             or die "define ${section}::chrLenFile first";
  my $chrFiles               = $config->{$section}{chrFiles}               or die "define ${section}::chrFiles first, it is the path of directory contains chromosome fasta files";
  my $ploidy                 = $config->{$section}{ploidy}                 or die "define ${section}::ploidy first, such as 2";
  my $coefficientOfVariation = $config->{$section}{coefficientOfVariation} or die "define ${section}::coefficientOfVariation first, such as 0.05";
  my $inputFormat            = $config->{$section}{inputFormat}            or die "define ${section}::inputFormat first";

  my $mateOrientation = $config->{$section}{mateOrientation}
    or die "define ${section}::mateOrientation first";    #0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)

  my $bedfile = parse_param_file( $config, $section, "bedfile", 0 );

  my $rawFiles = get_raw_files( $config, $section );
  my $groups;
  if(defined $config->{$section}{"groups"}){
    $groups = get_raw_files( $config, $section, "groups" );
  }else{
    $groups = {};
    for my $sampleName (sort keys %{$rawFiles}){
      $groups->{$sampleName} = $rawFiles->{$sampleName};
    } 
  }

  my $shfile = $pbsDir . "/${task_name}.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $groupName ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$groupName} };
    my $bamFile1;
    my $bamFile2;
    if ( scalar(@samples) > 1 ) {
      $bamFile2 = $rawFiles->{ $samples[0] }[0];    #control
      $bamFile1 = $rawFiles->{ $samples[1] }[0];    #sample
    }
    else {
      $bamFile1 = $rawFiles->{ $samples[0] }[0];    #sample
    }

    my $pbsFile = $self->pbsfile( $pbsDir, $groupName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $groupName );

    print SH "\$MYCMD ./$pbsName \n";
    my $curDir     = create_directory_or_die( $resultDir . "/$groupName" );
    my $configName = "${groupName}.conf";
    my $configFile = ${curDir} . "/$configName";

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

    open( CON, ">$configFile" ) or die $!;

    print CON "[general] 
chrLenFile=$chrLenFile 
ploidy=$ploidy 
coefficientOfVariation=$coefficientOfVariation 
chrFiles=$chrFiles
BedGraphOutput=TRUE

[sample]
mateFile=$bamFile1 
inputFormat=$inputFormat 
mateOrientation=$mateOrientation 

";
    if ( defined $bamFile2 ) {
      print CON "[control] 
mateFile=$bamFile2 
inputFormat=$inputFormat 
mateOrientation=$mateOrientation 

";
    }

    if ( defined $bedfile ) {
      print CON "[target]
captureRegions = $bedfile

";
    }

    close(CON);
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all freec tasks.\n";
}

sub performPileup {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $chrLenFile             = $config->{$section}{chrLenFile}             or die "define ${section}::chrLenFile first";
  my $chrFiles               = $config->{$section}{chrFiles}               or die "define ${section}::chrFiles first, it is the path of directory contains chromosome fasta files";
  my $ploidy                 = $config->{$section}{ploidy}                 or die "define ${section}::ploidy first, such as 2";
  my $coefficientOfVariation = $config->{$section}{coefficientOfVariation} or die "define ${section}::coefficientOfVariation first, such as 0.05";
  my $inputFormat            = $config->{$section}{inputFormat}            or die "define ${section}::inputFormat first";

  my $mateOrientation = $config->{$section}{mateOrientation}
    or die "define ${section}::mateOrientation first";    #0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs)

  my $snpFile = parse_param_file( $config, $section, "SNPfile", 1 );
  my $bedfile = parse_param_file( $config, $section, "bedfile", 0 );

  my $fastaFile = parse_param_file( $config, $section, "fasta_file", 1 );

  my $rawFiles = get_raw_files( $config, $section );
  my $groups = get_raw_files( $config, $section, "groups" );

  my $shfile = $pbsDir . "/${task_name}.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $groupName ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$groupName} };
    my $bamFile1;
    my $bamFile2;
    my $pileup1;
    my $pileup2;
    if ( scalar(@samples) > 1 ) {
      $bamFile2 = $rawFiles->{ $samples[0] }[0];    #control
      $pileup2  = $samples[0] . ".pileup.gz";
      $bamFile1 = $rawFiles->{ $samples[1] }[0];    #sample
      $pileup1  = $samples[1] . ".pileup.gz";
    }
    else {
      $bamFile1 = $rawFiles->{ $samples[0] }[0];    #sample
      $pileup1  = $samples[0] . ".pileup.gz";
    }

    my $pbsFile = $self->pbsfile( $pbsDir, $groupName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $groupName );

    print SH "\$MYCMD ./$pbsName \n";
    my $curDir     = create_directory_or_die( $resultDir . "/$groupName" );
    my $configName = "${groupName}.conf";
    my $configFile = ${curDir} . "/$configName";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir

echo freec_start=`date`
";

    print OUT "
if [ ! -s $pileup1 ]; then
  samtools mpileup -A -f $fastaFile $bamFile1 | gzip > $pileup1
fi
";

    if ( defined $pileup2 ) {
      print OUT "
if [ ! -s $pileup2 ]; then
  samtools mpileup -A -f $fastaFile $bamFile2 | gzip > $pileup2
fi
";
    }

    print OUT "
freec -conf $configName 

echo freec_end=`date`
";
    close OUT;
    print "$pbsFile created\n";

    open( CON, ">$configFile" ) or die $!;

    print CON "[general] 
chrLenFile=$chrLenFile 
ploidy=$ploidy 
coefficientOfVariation=$coefficientOfVariation 
chrFiles=$chrFiles
BedGraphOutput=TRUE

[sample]
mateFile=$pileup1 
inputFormat=pileup 
mateOrientation=$mateOrientation 

";
    if ( defined $bamFile2 ) {
      print CON "[control] 
mateFile=$pileup2 
inputFormat=pileup 
mateOrientation=$mateOrientation 

";
    }

    if ( defined $bedfile ) {
      print CON "[target]
captureRegions = $bedfile

";
    }

    print CON "[BAF]
SNPfile = $snpFile
";
    close(CON);
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
