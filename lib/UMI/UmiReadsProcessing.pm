#!/usr/bin/perl
package UMI::UmiReadsProcessing;

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
  $self->{_suffix} = "_umi";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );


  my $bwa_index = $config->{$section}{bwa_index};
  if ( !defined $bwa_index ) {
    $bwa_index = $config->{$section}{fasta_file} or die "define ${section}::bwa_index first";
  }
  #my $picard_jar = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1, not $self->using_docker());
  my $fgbio_jar = get_param_file( $config->{$section}{fgbio_jar}, "fgbio_jar", 1, not $self->using_docker());
  my $alignConsensusReads = get_option( $config, $section, "alignConsensusReads", 0 );
  #my $cutadapt = get_option( $config, $section, "cutadapt", 1 );
  
  my %mapped_bam_files = %{ get_raw_files( $config, $section ) };
  my %unmapped_bam_files = %{ get_raw_files( $config, $section,"unmapped_bam_files" ) };
  
  my $extension = get_option( $config, $section, "extension", "_umi" );
  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option || $java_option eq "" ) {
    $java_option = "-Xmx${memory}";
  }
  
  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";
  
  for my $sample_name ( sort keys %mapped_bam_files ) {
    my @sample_files = @{ $mapped_bam_files{$sample_name} };
    my $bam_file     = $sample_files[0];
    my $unmapped_bam_file = @{ $unmapped_bam_files{$sample_name} }[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );
    my $final_file;
    my $rmlist = "";

    print $pbs "echo Processing $sample_name \n"; 
    
# ##################################
# # align unmapped bam
# ##################################
#     print $pbs "echo 1. Align unmapped bam \n"; 
#     my $unmapped_fastq_prefix=basename($bam_file);
#     my $unmapped_fastq_file="${unmapped_fastq_prefix}.fastq.gz";
#     print $pbs "
# if [[ ! -s ${unmapped_fastq_file} ]]; then
#   java $java_option -jar $picard_jar SamToFastq I=$bam_file F=${unmapped_fastq_file} INTERLEAVE=true --COMPRESS_OUTPUTS_PER_RG
# fi
# ";
#     #cutadapt on fastq
#     if ($cutadapt) {
#       my $trimmed_fastq_file="${unmapped_fastq_prefix}.trimmed.fastq.gz";
#           print $pbs "
# if [[ ! -s ${trimmed_fastq_file} ]]; then
#   cutadapt --interleaved -n 2 -O 1 -q 20 -a AGATCGGAAGAGCACACGTC -A AGATCGGAAGAGCGTCGTGT -m 30 --trim-n -o ${trimmed_fastq_file} ${unmapped_fastq_file}
# fi
# ";
#       #change unmapped_fastq_prefix to trimmed file for next step
#       $unmapped_fastq_file=$trimmed_fastq_file;
#     }

#     my $mapped_sam_prefix=$sample_name.$extension;
#     print $pbs "
# if [[ ! -s ${mapped_sam_prefix}.sam ]]; then
#   bwa mem -M -p -t 8 $bwa_index ${unmapped_fastq_file} > ${mapped_sam_prefix}.sam
# fi
# ";

##################################
# merge alignment result with UMI tag
##################################
    my $umi_mapped_bam_prefix=basename($sample_name);
    my $umi_mapped_bam_merged_prefix=$umi_mapped_bam_prefix."_mappedAndMerged";
    $rmlist = $rmlist . " " . "${umi_mapped_bam_merged_prefix}.bam". " " . "${umi_mapped_bam_merged_prefix}.bai";
    print $pbs "
if [[ ! -s ${umi_mapped_bam_merged_prefix}.bam ]]; then
  echo 1. MergeBamAlignment 
  gatk MergeBamAlignment --UNMAPPED_BAM ${unmapped_bam_file} --ALIGNED_BAM ${bam_file} --OUTPUT ${umi_mapped_bam_merged_prefix}.bam -R $bwa_index \\
  -SO coordinate --ALIGNER_PROPER_PAIR_FLAGS true --ALIGNED_READS_ONLY true -MAX_GAPS -1 \\
  -ORIENTATIONS FR --VALIDATION_STRINGENCY SILENT --CREATE_INDEX true
fi
";

##################################
# generate consensus reads
##################################
    print $pbs "echo 2. Generate consensus reads \n"; 
    my $umi_mapped_grouped_bam_prefix=$umi_mapped_bam_merged_prefix."_grouped";
    $rmlist = $rmlist . " " . "${umi_mapped_grouped_bam_prefix}.bam";
    print $pbs "
if [[ ! -s ${umi_mapped_grouped_bam_prefix}.bam ]]; then
  java $java_option -jar $fgbio_jar GroupReadsByUmi --input=${umi_mapped_bam_merged_prefix}.bam --output=${umi_mapped_grouped_bam_prefix}.bam --strategy=adjacency --edits=1 --min-map-q=20 
fi
";
    my $umi_consensus_bam_prefix=$sample_name.$extension."_UnmappedConsensusReads";
    $rmlist = $rmlist . " " . "${umi_consensus_bam_prefix}.bam";
    print $pbs "
if [[ ! -s ${umi_consensus_bam_prefix}.bam ]]; then
  java $java_option -jar $fgbio_jar CallMolecularConsensusReads --input=${umi_mapped_grouped_bam_prefix}.bam --output=${umi_consensus_bam_prefix}.bam --error-rate-post-umi=30 --min-reads=1
fi
";
    my $umi_consensus_filtered_bam_prefix=$umi_consensus_bam_prefix."_PostFilter";
    $rmlist = $rmlist . " " . "${umi_consensus_filtered_bam_prefix}.bam";
    $final_file="${umi_consensus_filtered_bam_prefix}.sortByQuery.bam";
    print $pbs "
if [[ ! -s ${umi_consensus_filtered_bam_prefix}.bam ]]; then
  java $java_option -jar $fgbio_jar FilterConsensusReads --input=${umi_consensus_bam_prefix}.bam --output=${umi_consensus_filtered_bam_prefix}.bam --ref=$bwa_index --reverse-per-base-tags=true --min-reads=1 -E 0.05 -N 28 -e 0.1 -n 0.3
fi

#sort by query name so that SamToFastq won't have out of memory error
if [[ ! -s ${umi_consensus_filtered_bam_prefix}.sortByQuery.bam ]]; then
  sambamba sort -t 8 --sort-picard -o ${umi_consensus_filtered_bam_prefix}.sortByQuery.bam ${umi_consensus_filtered_bam_prefix}.bam
fi
";


##################################
# remap the filtered reads
##################################
if ($alignConsensusReads) {
    print $pbs "echo 3. Remap the filtered reads \n"; 
    my $unmapped_consensus_fastq_prefix=$sample_name.$extension."_consensus";
    $rmlist = $rmlist . " " . "${umi_consensus_filtered_bam_prefix}.sortByQuery.bam". " " . "${unmapped_consensus_fastq_prefix}.fastq";
    print $pbs "
if [[ ! -s ${unmapped_consensus_fastq_prefix}.fastq ]]; then
  gatk $java_option SamToFastq I=${umi_consensus_filtered_bam_prefix}.bam F=${unmapped_consensus_fastq_prefix}.fastq INTERLEAVE=true
fi
";
    my $mapped_consensus_sam_prefix=$unmapped_consensus_fastq_prefix;
    $rmlist = $rmlist . " " . "${mapped_consensus_sam_prefix}.sam";
    print $pbs "
if [[ ! -s ${mapped_consensus_sam_prefix}.sam ]]; then
  bwa mem -p -t 8 $bwa_index ${unmapped_consensus_fastq_prefix}.fastq > ${mapped_consensus_sam_prefix}.sam
fi
";
    my $umi_mapped_consensus_bam_prefix=$mapped_consensus_sam_prefix."_mapped";
     $final_file="${umi_mapped_consensus_bam_prefix}.bam";
    print $pbs "
if [[ ! -s ${umi_mapped_consensus_bam_prefix}.bam ]]; then
  gatk $java_option MergeBamAlignment UNMAPPED=${umi_consensus_filtered_bam_prefix}.bam ALIGNED=${mapped_consensus_sam_prefix}.sam O=${umi_mapped_consensus_bam_prefix}.bam R=$bwa_index SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true 
fi
";
}

  if ($rmlist ne "") {
              print $pbs "
if [ -s $final_file ]; then
  rm $rmlist
fi
";
  } 
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }
  
  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );
  my $result = {};

  my $extension = get_option( $config, $section, "extension", "_umi" );
  my $alignConsensusReads = get_option( $config, $section, "alignConsensusReads", 1 );

  my %unmapped_bam_files = %{ get_raw_files( $config, $section ) };
  for my $sample_name ( sort keys %unmapped_bam_files ) {
    my $bamOut="";
    if ($alignConsensusReads) {
      $bamOut       = $sample_name.$extension."_consensus_mapped.bam";
    } else {
      $bamOut       = $sample_name.$extension."_UnmappedConsensusReads_PostFilter.sortByQuery.bam";
    }
    
    my @result_files = ();
    push( @result_files, "${result_dir}/${bamOut}" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
