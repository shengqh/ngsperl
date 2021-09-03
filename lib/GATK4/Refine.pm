#!/usr/bin/perl
package GATK4::Refine;

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
  $self->{_use_tmp_folder} = 1;
  $self->{_docker_prefix}   = "gatk4_";
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
    my $vcfFiles     = $config->{$section}{indel_vcf_files};
    my @vcfFileArray = is_array($vcfFiles) ? @{$vcfFiles} : ($vcfFiles);
    my @out          = keys %{ { map { ( $_ => 1 ) } ( @sitesVcfFiles, @vcfFileArray ) } };
    @sitesVcfFiles = @out;
  }

  my $result = "";
  foreach my $vcf (@sitesVcfFiles) {
    $result = $result . " -knownSites $vcf";
  }

  return $result;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my $knownsitesvcf = getKnownSitesVcf( $config, $section );

  my $java_option = $self->get_java_option( $config, $section, $memory );
  
  #other options
  my $faFile     = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $fixMisencodedQuals = get_option( $config, $section, "fixMisencodedQuals", 0 ) ? "--fixMisencodedQuals" : "";
  my $sorted           = get_option( $config, $section, "sorted", 0 );
  my $remove_duplicate = get_option( $config, $section, "remove_duplicate",         1 );
  my $mark_duplicate   = get_option( $config, $section, "mark_duplicate",           0 );

  if(! $mark_duplicate) {
    $remove_duplicate = 1;
  }

  my $removeDupLabel = $remove_duplicate ? "true" : "false";
  my $rmdupResultName = $remove_duplicate ? ".rmdup" : ".markdup";

  my $indel_vcf = "";
  my $vcfFiles = $config->{$section}{indel_vcf_files} or die "Define indel_vcf_files in section $section first.";
  my @vcfFileArray = is_array($vcfFiles) ? @{$vcfFiles} : ($vcfFiles);
  foreach my $vcf (@vcfFileArray) {
    $indel_vcf = $indel_vcf . " --known-sites $vcf";
  }

  my $restrict_intervals = "";
  my $bedFile = get_param_file( $config->{$section}{bed_file}, "bed_file", 0 );
  if ( defined $bedFile and $bedFile ne "" ) {
    $restrict_intervals = "-L $bedFile";
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sampleFile   = $sample_files[0];

    my $presortedFile = "${sample_name}.sorted.bam";

    my $rmdupFile = "${sample_name}${rmdupResultName}.bam";
    my $rmdupFileIndex = change_extension( $rmdupFile, ".bai" );

    my $recalibration_report_filename = "${sample_name}${rmdupResultName}.recal_data.csv";

    my $final_file = "${sample_name}${rmdupResultName}.recal.bam";
    my $final_file_md5 = "${sample_name}${rmdupResultName}.recal.bam.md5";
    my $final_file_index = change_extension( $final_file, ".bai" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "if [[ ! -s $result_dir/$final_file ]]; then
  \$MYCMD ./$pbs_name 
fi
";

    my $log_desc = $cluster->get_log_description($log);

    my $rmlist       = "";
    my $resultPrefix = $sample_name;
    my $inputFile    = $sample_files[0];

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command, 0, $inputFile );

    my $localized_files = [];
    @sample_files = @{$self->localize_files_in_tmp_folder($pbs, \@sample_files, $localized_files, [".bai"])};
    $inputFile    = $sample_files[0];
    my $inputFileIndex    = "${inputFile}.bai";

    if ( !$sorted ) {
      print $pbs "
if [[ (-s $inputFile) && ! (-s $presortedFile) ]]; then
  echo sort=`date` 
  samtools sort -m 4G -o $presortedFile $inputFile
  samtools index $presortedFile
fi
";
      $inputFile = $presortedFile;
      $inputFileIndex = "${presortedFile}.bai";
      $rmlist = "$rmlist $presortedFile ${presortedFile}.bai";
    }

    my $lastFileList = "";

    print $pbs "

echo MarkDuplicates=`date` 
gatk --java-options \"$java_option\" \\
  MarkDuplicates \\
  --INPUT $inputFile \\
  --OUTPUT $rmdupFile \\
  --METRICS_FILE ${rmdupFile}.metrics \\
  --VALIDATION_STRINGENCY SILENT \\
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
  --REMOVE_DUPLICATES $removeDupLabel \\
  --CREATE_INDEX true \\
  --ASSUME_SORTED true \\
  --CREATE_MD5_FILE false

status=\$?
if [[ \$status -ne 0 ]]; then
  touch $sample_name.MarkDuplicates.failed
  rm -f $rmdupFile ${rmdupFile}.metrics $rmdupFileIndex
else
  echo BaseRecalibrator=`date` 
  gatk --java-options \"$java_option\" \\
    BaseRecalibrator \\
    -R $faFile \\
    -I $rmdupFile \\
    --use-original-qualities $fixMisencodedQuals \\
    $indel_vcf \\
    -O $recalibration_report_filename $restrict_intervals 

  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sample_name.BaseRecalibrator.failed
  else
    echo ApplyBQSR=`date`
    gatk --java-options \"$java_option\" \\
      ApplyBQSR \\
      -R $faFile \\
      -I $rmdupFile \\
      -O $final_file \\
      -bqsr $recalibration_report_filename \\
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \\
      --add-output-sam-program-record \\
      --create-output-bam-md5 \\
      --use-original-qualities

    status=\$?
    if [[ \$status -ne 0 ]]; then
      touch $sample_name.ApplyBQSR.failed
    else
      touch $sample_name.succeed
      ln $final_file_index ${final_file}.bai

      echo flagstat = `date` 
      samtools flagstat $final_file > ${final_file}.stat 

      echo flagstat = `date` 
      samtools idxstats $final_file > ${final_file}.chromosome.count 

      rm $result_dir/$sample_name.*.failed
      rm $rmlist $rmdupFile $rmdupFileIndex
    fi
  fi
fi

";

    $self->clean_temp_files($pbs, $localized_files);

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all GATK refine tasks. \n ";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $remove_duplicate = get_option( $config, $section, "remove_duplicate",         1 );
  my $mark_duplicate   = get_option( $config, $section, "mark_duplicate",           0 );

  if(! $mark_duplicate) {
    $remove_duplicate = 1;
  }

  my $rmdupResultName = $remove_duplicate ? ".rmdup" : ".markdup";

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $final_file = "${sample_name}${rmdupResultName}.recal.bam";
    my $final_file_index = $final_file . ".bai";

    my @result_files = ();
    push( @result_files, "${result_dir}/${final_file}" );
    push( @result_files, "${result_dir}/${final_file_index}" );
    push( @result_files, "${result_dir}/${final_file}.chromosome.count" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
