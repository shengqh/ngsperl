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

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(cnvnator conifer cnmops freec)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub cnvnator {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );
  my $binsize        = $config->{$section}{binsize}        or die "define ${section}::binsize first";
  my $chromosome_dir = $config->{$section}{chromosome_dir} or die "define ${section}::chromosome_dir first";

  my $genome    = $config->{$section}{genome};
  my $genomestr = "";

  if ( defined $genome ) {
    $genomestr = "-genome " . $genome;
  }

  my $isbamsorted = $config->{$section}{isbamsorted};
  if ( !defined($isbamsorted) ) {
    $isbamsorted = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $bamFile = $sampleFiles[0];

    if ( !$isbamsorted ) {
      ( $bamFile, my $bamSorted ) = get_sorted_bam($bamFile);
    }
    my $pbsName = "cnvnator_${sampleName}.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/cnvnator_${sampleName}.log";
    my $curDir   = create_directory_or_die( $resultDir . "/$sampleName" );
    my $rootFile = $sampleName . ".root";
    my $callFile = $sampleName . ".call";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

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

  #`qsub $pbsFile`;
}

sub conifer {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );
  my $conifer = $config->{$section}{conifer} or die "define conifer program location first.\nconifer => \"location\"";
  my $bedfile = $config->{$section}{bedfile} or die "define $section:bedfile first";

  my $isbamsorted = $config->{$section}{isbamsorted};
  if ( !defined($isbamsorted) ) {
    $isbamsorted = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}_rpkm.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  create_directory_or_die( $resultDir . "/rpkm" );
  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $pbsName = "${sampleName}_rpkm.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${sampleName}_rpkm.log";
    my $bamFile = $sampleFiles[0];

    if ( !$isbamsorted ) {
      ( $bamFile, my $bamSorted ) = get_sorted_bam($bamFile);

      #print $bamFile . "\n";
    }
    my $rpkm = "rpkm/" . $sampleName . ".rpkm";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir

echo rpkm=`date`

if [ ! -s $rpkm ]; then
  echo conifer=`date`
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

  my $pbsName   = "${task_name}_after_rpkm.pbs";
  my $pbsFile   = "${pbsDir}/$pbsName";
  my $hdf5File  = "${task_name}.hdf5";
  my $svalsFile = "${task_name}.svals";
  my $plotFile  = "${task_name}.png";
  my $sdFile  = "${task_name}.sd";
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

  #`qsub $pbsFile`;
}

sub cnmops {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );
  my $bedfile = $config->{$section}{bedfile};

  my $isbamsorted = $config->{$section}{isbamsorted};
  if ( !defined($isbamsorted) ) {
    $isbamsorted = 0;
  }

  my $pairmode = $config->{$section}{pairmode};
  if ( !defined($pairmode) ) {
    $pairmode = "unpaired";
  }

  my $callFile = "${task_name}.call";

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $rfile = $resultDir . "/cnmops_${task_name}.r";
  open( R, ">$rfile" ) or die "Cannot create $rfile";
  print R "library(cn.mops) \n";
  print R "setwd(\"$resultDir\") \n";
  print R "SampleNames <- c( \n";
  my $isfirst = 1;
  for my $sampleName ( sort keys %rawFiles ) {
    if ($isfirst) {
      print R "\"$sampleName\"\n";
      $isfirst = 0;
    }
    else {
      print R ",\"$sampleName\"\n";
    }
  }
  print R ") \n";
  print R "BAMFiles <- c( \n";

  $isfirst = 1;
  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $bamFile     = $sampleFiles[0];

    if ( !$isbamsorted ) {
      ( $bamFile, my $bamSorted ) = get_sorted_bam($bamFile);
    }

    if ($isfirst) {
      print R "\"$bamFile\"\n";
      $isfirst = 0;
    }
    else {
      print R ",\"$bamFile\"\n";
    }
  }
  print R ") \n";

  if ( defined $bedfile ) {
    print R "
segments <- read.table(\"$bedfile\", sep=\"\\t\", as.is=TRUE, header=T) 
gr <- GRanges(segments[,1], IRanges(segments[,2],segments[,3]), gene=segments[,4]) 
x <- getSegmentReadCountsFromBAM(BAMFiles, GR=gr, sampleNames=SampleNames, mode=\"$pairmode\") 
save(x, file=\"${task_name}_x_getSegmentReadCountsFromBAM.Rdata\") 
resCNMOPS <- exomecn.mops(x) 
save(resCNMOPS, file=\"${task_name}_resCNMOPS_exomecn.mops.Rdata\")
";
  }
  else {
    print R "
x <- getReadCountsFromBAM(BAMFiles, sampleNames=SampleNames, mode=\"$pairmode\") 
save(x, file=\"${task_name}_x_getReadCountsFromBAM.Rdata\") 
resCNMOPS <- cn.mops(x) 
save(resCNMOPS, file=\"${task_name}_resCNMOPS_cn.mops.Rdata\") 
";
  }

  print R "
cnvs<-resCNMOPS\@cnvs 
d<-cbind(substring(as.character(cnvs\@elementMetadata\@listData\$sampleName),2), 
         as.character(cnvs\@seqnames),  
         as.character(cnvs\@ranges\@start), 
         as.character(as.numeric(cnvs\@ranges\@start) + as.numeric(cnvs\@ranges\@width) - 1), 
         as.character(cnvs\@ranges\@width), 
         as.character(cnvs\@elementMetadata\@listData\$CN), 
         as.character(cnvs\@elementMetadata\@listData\$median), 
         as.character(cnvs\@elementMetadata\@listData\$mean)) 
colnames(d)<-c(\"sample\",\"chr\",\"start\",\"end\", \"length\",\"type\",\"median\",\"mean\") 
d<-d[order(d[,\"sample\"], as.numeric(d[,\"chr\"]), as.numeric(d[,\"start\"])),] 
d[,\"chr\"]<-paste0(\"chr\",d[,\"chr\"]) 
d[,\"type\"]<-apply(d,1,function(x){ 
  if(as.numeric(x[\"median\"]) < 0){ 
    return (\"DELETION\") 
  }else{ 
    return (\"DUPLICATION\") 
  } 
}) 
write.table(d, file=\"$callFile\",sep=\"\\t\",col.names=T,row.names=F,quote=F) 
";
  close R;

  my $pbsFile = "${pbsDir}/cnmops_${task_name}.pbs";
  my $log     = "${logDir}/cnmops_${task_name}.log";

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $pbsDir
echo cnmops=`date`
R --vanilla < $rfile 
echo finished=`date`
";
  close OUT;

  print "$pbsFile created\n";
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
