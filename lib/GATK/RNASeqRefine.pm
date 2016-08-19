#!/usr/bin/perl
package GATK::RNASeqRefine;

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
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_rrf";
  bless $self, $class;
  return $self;
}

#RNASeq: based on https://software.broadinstitute.org/gatk/guide/article?id=3891
sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my @vcfFiles = @{ $config->{$section}{vcf_files} } or die "Define vcf_files in section $section first.";
  my $gatk_jar   = get_param_file( $config->{$section}{gatk_jar},   "gatk_jar",   1 );
  my $picard_jar = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1 );

  my $gatk_option = get_option( $config, $section, "gatk_option", "" );

  my $sorted = get_option( $config, $section, "sorted", 0 );

  my $replaceReadGroup   = get_option( $config, $section, "replace_read_group",       0 );
  my $reorderChromosome  = get_option( $config, $section, "reorder_chromosome",       0 );
  my $fixMisencodedQuals = get_option( $config, $section, "fixMisencodedQuals",       0 ) ? "-fixMisencodedQuals" : "";
  my $slimPrintReads     = get_option( $config, $section, "slim_print_reads",         1 );
  my $baq                = get_option( $config, $section, "samtools_baq_calibration", 0 );

  my $knownvcf      = "";
  my $knownsitesvcf = "";

  foreach my $vcf (@vcfFiles) {
    $knownvcf      = $knownvcf . " -known $vcf";
    $knownsitesvcf = $knownsitesvcf . " -knownSites $vcf";
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files   = @{ $raw_files{$sample_name} };
    my $sampleFile     = $sample_files[0];
    my $sampleFileName = basename($sampleFile);

    my $rmFiles       = "";
    my $inputFile     = $sampleFile;
    my $presortedFile = "";
    my $sortCmd       = "";
    if ( !$sorted ) {
      my $presortedPrefix = $sample_name . ".sorted";
      $presortedFile = $presortedPrefix . ".bam";
      $sortCmd       = "samtools sort -@ $thread -m 4G $sampleFile $presortedPrefix";
      $inputFile     = $presortedFile;
      $rmFiles       = $presortedFile;
    }

    my $rmdupFile = $sample_name . ".rmdup.bam";

    my $splitInput   = $rmdupFile;
    my $replaceCmd   = "";
    my $replacedFile = "";
    if ($replaceReadGroup) {
      $replacedFile = $sample_name . ".rmdup.rgreplaced.bam";
      $replaceCmd =
"if [ ! -s $replacedFile ]; then java -jar $picard_jar AddOrReplaceReadGroups I=$rmdupFile O=$replacedFile ID=$sample_name LB=$sample_name SM=$sample_name PL=ILLUMINA PU=ILLUMINA; samtools index $replacedFile; fi;";
      $splitInput = $replacedFile;
      $rmFiles    = $rmFiles . " " . $replacedFile . " " . $replacedFile . ".bai";
    }

    my $reorderCmd  = "";
    my $reorderFile = "";
    if ($reorderChromosome) {
      $reorderFile = $sample_name . ".rmdup.reorder.bam";
      $reorderCmd  = "if [ ! -s $reorderFile ]; then java -jar $picard_jar ReorderSam I=$splitInput O=$reorderFile REFERENCE=$faFile; samtools index $reorderFile; fi;";
      $splitInput  = $reorderFile;
      $rmFiles     = $rmFiles . " " . $reorderFile . " " . $reorderFile . ".bai";
    }

    my $splitFile = $sample_name . ".rmdup.split.bam";
    my $grpFile   = $splitFile . ".grp";
    my $recalFile = $sample_name . ".rmdup.split.recal.bam";

    my $final_file = $recalFile;

    my $rmlist = "";

    my $slimcmd = "";
    if ($slimPrintReads) {
      $rmlist = "$final_file ${final_file}.bai";
      my $slimFile = $sample_name . ".rmdup.split.recal.slim.bam";
      $slimcmd = "if [[ -s $final_file && ! -s $slimFile ]]; then
  echo slim=`date` 
  samtools view -h $final_file | sed 's/\\tBD\:Z\:[^\\t]*//' | sed 's/\\tPG\:Z\:[^\\t]*//' | sed 's/\\tBI\:Z\:[^\\t]*//' | samtools view -S -b > $slimFile
  samtools index $slimFile
fi
";
      my $final_file = $slimFile;
    }

    my $baqcmd = "";
    if ($baq) {
      $rmlist = $rmlist . " $final_file ${final_file}.bai";
      my $baq_file = $sample_name . ".rmdup.split.recal.slim.baq.bam";
      $baqcmd = "
if [[ -s $final_file && ! -s $baq_file ]]; then
  echo baq=`date` 
  samtools calmd -Abr $final_file $faFile > $baq_file
  samtools index $baq_file
fi      
";
      $final_file = $baq_file;
    }

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    print $pbs "
if [ ! -s $rmdupFile ]; then
  echo MarkDuplicates=`date` 
  $sortCmd
  java $option -jar $picard_jar MarkDuplicates I=$inputFile O=$rmdupFile ASSUME_SORTED=true REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${rmdupFile}.metrics
fi

if [ ! -s $splitFile ]; then
  echo SplitNCigarReads=`date` 
  $replaceCmd
  $reorderCmd
  java $option -jar $gatk_jar -T SplitNCigarReads -R $faFile -I $splitInput -o $splitFile -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS $fixMisencodedQuals
fi

if [[ -s $splitFile && ! -s $grpFile ]]; then
  echo BaseRecalibrator=`date` 
  java $option -jar $gatk_jar -T BaseRecalibrator -nct $thread -rf BadCigar -R $faFile -I $splitFile $knownsitesvcf -o $grpFile
fi

if [[ -s $splitFile && -s $grpFile && ! -s $recalFile ]]; then
  echo PrintReads=`date`
  java $option -jar $gatk_jar -T PrintReads -nct $thread -rf BadCigar -R $faFile -I $splitFile -BQSR $grpFile -o $recalFile 
fi

$slimcmd

$baqcmd

if [[ -s $final_file && ! -s ${final_file}.stat ]]; then
  echo flagstat=`date` 
  samtools flagstat $final_file > ${final_file}.stat
  rm $rmFiles $rmdupFile ${sample_name}.rmdup.bai ${rmdupFile}.metrics $splitFile ${sample_name}.rmdup.split.bai $grpFile $rmlist
fi
  
";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all GATK refine tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $slimPrintReads = get_option( $config, $section, "slim_print_reads",         1 );
  my $baq            = get_option( $config, $section, "samtools_baq_calibration", 0 );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $recalFile  = $sample_name . ".rmdup.split.recal.bam";
    my $final_file = $recalFile;
    if ($slimPrintReads) {
      my $slimFile   = $sample_name . ".rmdup.split.recal.slim.bam";
      my $final_file = $slimFile;
    }

    if ($baq) {
      my $baq_file = $sample_name . ".rmdup.split.recal.slim.baq.bam";
      $final_file = $baq_file;
    }
    my @result_files = ();
    push( @result_files, "${result_dir}/${sample_name}/${final_file}" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
