#!/usr/bin/perl
package CQS::GATKRefine;

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
  $self->{_name} = "GATKRefine";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my @vcfFiles = @{ $config->{$section}{vcf_files} } or die "Define vcf_files in section $section first.";
  my $gatk_jar           = get_param_file( $config->{$section}{gatk_jar},           "gatk_jar",           1 );
  my $markDuplicates_jar = get_param_file( $config->{$section}{markDuplicates_jar}, "markDuplicates_jar", 1 );
  my $thread_count       = $config->{$section}{thread_count};

  my $knownvcf      = "";
  my $knownsitesvcf = "";

  foreach my $vcf (@vcfFiles) {
    $knownvcf      = $knownvcf . " -known $vcf";
    $knownsitesvcf = $knownsitesvcf . " -knownSites $vcf";
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}_rf.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  if ($sh_direct) {
    print SH "export MYCMD=\"bash\" \n";
  }
  else {
    print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
  }

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles   = @{ $rawFiles{$sampleName} };
    my $sampleFile    = $sampleFiles[0];
    my $sampleFileName = basename($sampleFile);
    
    my $intervalFile  = $sampleFileName . ".intervals";
    my $realignedFile = change_extension( $sampleFileName, ".realigned.bam" );
    my $grpFile       = $realignedFile . ".grp";
    my $recalFile     = change_extension( $realignedFile, ".recal.bam" );
    my $rmdupFile     = change_extension( $recalFile, ".rmdup.bam" );
    my $sortedPrefix  = change_extension( $rmdupFile, "_sorted" );
    my $sortedFile    = $sortedPrefix . ".bam";

    my $pbsName = "${sampleName}_rf.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";
    my $curDir  = create_directory_or_die( $resultDir . "/$sampleName" );
    my $log     = "${logDir}/${sampleName}_rf.log";

    print SH "\$MYCMD ./$pbsName \n";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $curDir

if [ ! -s $intervalFile ]; then
  echo RealignerTargetCreator=`date` 
  java $option -jar $gatk_jar -T RealignerTargetCreator -I $sampleFile -R $faFile $knownvcf -nt $thread_count -o $intervalFile
fi

if [[ -s $intervalFile && ! -s $realignedFile ]]; then
  echo IndelRealigner=`date` 
  #InDel parameter referenced: http://www.broadinstitute.org/gatk/guide/tagged?tag=local%20realignment
  java $option -Djava.io.tmpdir=tmpdir -jar $gatk_jar -T IndelRealigner -I $sampleFile -R $faFile -targetIntervals $intervalFile $knownvcf --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4 -o $realignedFile 
fi

if [[ -s $realignedFile && ! -s $grpFile ]]; then
  echo BaseRecalibrator=`date` 
  java $option -jar $gatk_jar -T BaseRecalibrator -rf BadCigar -R $faFile -I $realignedFile $knownsitesvcf -o $grpFile -plots ${grpFile}.pdf
fi

if [[ -s $grpFile && ! -s $recalFile ]]; then
  echo PrintReads=`date`
  java $option -jar $gatk_jar -T PrintReads -rf BadCigar -R $faFile -I $realignedFile -BQSR $grpFile -o $recalFile 
fi

if [[ -s $recalFile && ! -s $rmdupFile ]]; then
  echo RemoveDuplicate=`date` 
  java $option -jar $markDuplicates_jar I=$recalFile O=$rmdupFile M=${rmdupFile}.matrix VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true
fi

if [[ -s $rmdupFile && ! -s $sortedFile ]]; then
  echo BamSort=`date` 
  samtools sort $rmdupFile $sortedPrefix 
fi

if [[ -s $sortedFile && ! -s ${sortedFile}.bai ]]; then
  echo BamIndex=`date` 
  samtools index $sortedFile
fi
  
echo finished=`date`

exit 1;
";

    close OUT;

    print "$pbsFile created\n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all GATK refine tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $sampleFile  = $sampleFiles[0];
    my $sampleFileName = basename($sampleFile);
    my $sortedFile  = change_extension( $sampleFileName, ".realigned.recal.rmdup_sorted.bam" );
    my @resultFiles = ();
    push( @resultFiles, "${resultDir}/${sampleName}/${sortedFile}" );
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
