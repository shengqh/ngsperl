#!/usr/bin/perl
package GATK::Refine;

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
  $self->{_suffix} = "_rf";
  bless $self, $class;
  return $self;
}

sub getKnownSitesVcf {
  my ( $config, $section ) = @_;

  my $sitesVcfFiles;
  if ( defined $config->{$section}{known_vcf_files} ) {
    $sitesVcfFiles = $config->{$section}{known_vcf_files};
  }
  else {
    $sitesVcfFiles = $config->{$section}{vcf_files} or die "Define vcf_files in section $section first.";
  }
  my @sitesVcfFiles = @{$sitesVcfFiles};

  if ( defined $config->{$section}{indel_vcf_files} ) {
    my $vcfFiles = $config->{$section}{indel_vcf_files};
    my @vcfFiles = @{$vcfFiles};
    my @out      = keys %{ { map { ( $_ => 1 ) } ( @sitesVcfFiles, @vcfFiles ) } };
    @sitesVcfFiles = @out;
  }

  my $result = "";
  foreach my $vcf (@sitesVcfFiles) {
    $result = $result . " -knownSites $vcf";
  }

  return $result;
}

#https://software.broadinstitute.org/gatk/blog?id=7712
#indel local alignment has been dropped out, since it is just not useful enough anymore when you're calling variants with haplotype-based tools like HaplotypeCaller and MuTect2
#
#https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_engine_CommandLineGATK.php--simplifyBAM
#--simplifyBAM / -simplifyBAM
#Strip down read content and tags
#If provided, output BAM/CRAM files will be simplified to include only key reads for downstream variation discovery analyses (removing duplicates, PF-, non-primary reads), as well stripping all extended tags from the kept reads except the read group identifier
#Based on this description, it should be safe to remove BD:Z, PG:Z and BI:Z in bam file.
sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = get_parameter( $config, $section );

  my $knownsitesvcf = getKnownSitesVcf( $config, $section );

  #other options
  my $faFile     = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $gatk_jar   = get_param_file( $config->{$section}{gatk_jar},   "gatk_jar",   1 );
  my $picard_jar = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1 );

  my $fixMisencodedQuals = get_option( $config, $section, "fixMisencodedQuals", 0 ) ? "-fixMisencodedQuals" : "";

  my $sorted               = get_option( $config, $section, "sorted",                   0 );
  my $remove_duplicate     = get_option( $config, $section, "remove_duplicate",         1 );
  my $indelRealignment     = get_option( $config, $section, "indel_realignment",        0 );
  my $slim                 = get_option( $config, $section, "slim_print_reads",         1 );
  my $use_self_slim_method = get_option( $config, $section, "use_self_slim_method",     0 );
  my $baq                  = get_option( $config, $section, "samtools_baq_calibration", 0 );

  my $indel_vcf = "";
  if ($indelRealignment) {
    my $vcfFiles = $config->{$section}{indel_vcf_files} or die "Define indel_vcf_files in section $section first.";
    foreach my $vcf (@$vcfFiles) {
      $indel_vcf = $indel_vcf . " -known $vcf";
    }
  }

  my $bedFile = get_param_file( $config->{$section}{bed_file}, "bed_file", 0 );
  my $restrict_intervals = "";

  if ( defined $bedFile and $bedFile ne "" ) {
    $restrict_intervals = "-L $bedFile";
    my $interval_padding = get_option( $config, $section, "interval_padding", 0 );
    if ( defined $interval_padding and $interval_padding != 0 ) {
      $restrict_intervals = $restrict_intervals . " -ip $interval_padding";
    }
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sampleFile   = $sample_files[0];

    my $rmdupResultName = $remove_duplicate ? ".rmdup" : "";
    my $slimResultName  = $slim             ? ".slim"  : "";
    my $indelResultName = $indelRealignment ? ".indel" : "";
    my $baqResultName   = $baq              ? ".baq"   : "";
    my $final_file      = "${sample_name}${rmdupResultName}.recal${indelResultName}${slimResultName}${baqResultName}.bam";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command );

    my $rmlist       = "";
    my $inputFile    = $sampleFile;
    my $resultPrefix = $sample_name;
    if ( !$sorted ) {
      my $presortedFile = $sample_name . ".sorted.bam";
      print $pbs "
if [[ -s $inputFile && ! -s $presortedFile ]]; then
  echo sort=`date` 
  samtools sort -@ $thread -m 4G -o $presortedFile $inputFile
  samtools index $presortedFile
fi
";
      $inputFile = $presortedFile;
      $rmlist    = $rmlist . " $presortedFile ${presortedFile}.bai";
    }

    if ($remove_duplicate) {
      my $rmdupFile = $sample_name . ".rmdup.bam";
      my $rmdupFileIndex = change_extension( $rmdupFile, ".bai" );
      print $pbs "
if [[ -s $inputFile && ! -s $rmdupFile ]]; then
  echo MarkDuplicates=`date` 
  java $option -jar $picard_jar MarkDuplicates I=$inputFile O=$rmdupFile ASSUME_SORTED=true REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${rmdupFile}.metrics
fi
";
      $inputFile = $rmdupFile;
      $rmlist    = $rmlist . " $rmdupFile $rmdupFileIndex";
    }

    my $recalTable   = "${sample_name}${rmdupResultName}.recal.table";
    my $recalFile    = "${sample_name}${rmdupResultName}.recal.bam";
    my $slimFile     = "${sample_name}${rmdupResultName}.recal${slimResultName}.bam";
    my $printOptions = "";
    if ( $slim and !$use_self_slim_method ) {
      $printOptions = " --simplifyBAM";
      $recalFile    = $slimFile;
    }
    my $recalFileIndex = change_extension( $recalFile, ".bai" );

    print $pbs "
if [[ -s $inputFile && ! -s $recalTable ]]; then
  echo BaseRecalibrator=`date` 
  java $option -jar $gatk_jar -T BaseRecalibrator -nct $thread -rf BadCigar -R $faFile -I $inputFile $knownsitesvcf -o $recalTable $restrict_intervals
fi

if [[ -s $recalTable && ! -s $recalFile ]]; then
  echo PrintReads=`date`
  java $option -jar $gatk_jar -T PrintReads $printOptions -nct $thread -rf BadCigar -R $faFile -I $inputFile -BQSR $recalTable -o $recalFile 
fi
";
    if ( $slim and $use_self_slim_method ) {
      $rmlist = $rmlist . " $recalFile $recalFile.bai";
      print $pbs "
if [[ -s $recalFile && ! -s $slimFile ]]; then
  echo slim=`date` 
  samtools view -h $recalFile | sed 's/\\tBD\:Z\:[^\\t]*//' | sed 's/\\tPG\:Z\:[^\\t]*//' | sed 's/\\tBI\:Z\:[^\\t]*//' | samtools view -S -b > $slimFile
  samtools index $slimFile
fi  
";
    }
    $inputFile = $slimFile;

    if ($indelRealignment) {
      my $intervalFile = "${sample_name}${rmdupResultName}.recal${slimResultName}${indelResultName}.intervals";
      my $indel_file   = "${sample_name}${rmdupResultName}.recal${slimResultName}${indelResultName}.bam";
      $rmlist = $rmlist . " $inputFile $inputFile.bai $intervalFile";
      print $pbs "
if [[ -s $inputFile && ! -s $indel_file ]]; then
  echo RealignerTargetCreator=`date` 
  java $option -jar $gatk_jar -T RealignerTargetCreator -nt $thread $fixMisencodedQuals -I $inputFile -R $faFile $indel_vcf -o $intervalFile $restrict_intervals
fi

if [[ -s $intervalFile && ! -s $indel_file ]]; then
  echo IndelRealigner=`date` 
  #InDel parameter referenced: http://www.broadinstitute.org/gatk/guide/tagged?tag=local%20realignment
  java $option -Djava.io.tmpdir=tmpdir -jar $gatk_jar -T IndelRealigner $fixMisencodedQuals -I $slimFile -R $faFile -targetIntervals $intervalFile $indel_vcf --consensusDeterminationModel USE_READS -LOD 0.4 -o $indel_file 
fi  
";
      $inputFile = $indel_file;
    }

    if ($baq) {
      print $pbs "
if [[ -s $inputFile && !-s $final_file ]]; then
  echo baq = `date` 
  samtools calmd -Abr $inputFile $faFile > $final_file samtools index $final_file 
fi
";
      $rmlist = $rmlist . " $inputFile $inputFile.bai ";
    }

    print $pbs "
if [[ -s $final_file && !-s ${final_file}.stat ]]; then 
  echo flagstat = `date` 
  samtools flagstat $final_file > ${final_file}.stat 
  rm $rmlist 
fi
";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all GATK refine tasks . \n ";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $remove_duplicate = get_option( $config, $section, "remove_duplicate",         1 );
  my $indelRealignment = get_option( $config, $section, "indel_realignment",        0 );
  my $slim             = get_option( $config, $section, "slim_print_reads",         1 );
  my $baq              = get_option( $config, $section, "samtools_baq_calibration", 0 );

  my $rmdupResultName = $remove_duplicate ? ".rmdup" : "";
  my $indelResultName = $indelRealignment ? ".indel" : "";
  my $slimResultName  = $slim             ? ".slim"  : "";
  my $baqResultName   = $baq              ? ".baq"   : "";

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $final_file   = "${sample_name}${rmdupResultName}.recal${slimResultName}${indelResultName}${baqResultName}.bam ";
    my @result_files = ();
    push( @result_files, "${result_dir}/${final_file}" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
