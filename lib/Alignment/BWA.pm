#!/usr/bin/perl
package Alignment::BWA;

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
use Alignment::AlignmentUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_bwa";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $sort_memory = $thread == 1 ? $memory : "4G";

  my $selfname = $self->{_name};

  my $cleansam                = get_option( $config, $section, "cleansam",                0 );
  my $chromosome_grep_pattern = get_option( $config, $section, "chromosome_grep_pattern", "" );
  my $sortByCoordinate        = get_option( $config, $section, "sort_by_coordinate",      1 );
  my $mark_duplicates         = hasMarkDuplicate( $config->{$section} );

  $option = $option . " -M";

  if ( !( $option =~ /\s-t\s/ ) ) {
    if ( $thread > 1 ) {
      $option = $option . " -t " . $thread;
    }
  }

  my $bwa_index = $config->{$section}{bwa_index};
  if ( !defined $bwa_index ) {
    $bwa_index = $config->{$section}{fasta_file} or die "define ${section}::bwa_index first";
  }

  my $picard_jar = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1, not $self->using_docker() );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sample_files_str = ( scalar(@sample_files) == 2 ) ? "\"" . $sample_files[0] . "\" \"" . $sample_files[1] . "\"" : "\"" . $sample_files[0] . "\"";

    my $unsorted_bam_file = $sample_name . ".unsorted.bam";
    my $clean_bam_file    = $sample_name . ".unsorted.clean.bam";
    my $rmdup_bam_file    = $sample_name . ".unsorted.clean.rmdup.bam";
    my $bam_file          = $sample_name . ".bam";
    my $tag               = get_bam_tag($sample_name);

    my $rg = "\@RG\\tID:${sample_name}\\tPU:illumina\\tLB:${sample_name}\\tSM:${sample_name}\\tPL:illumina";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $final_file = $sample_name . ( $mark_duplicates ? ".rmdup" : "" ) . ( $sortByCoordinate ? ".sorted" : ".unsorted" ) . ".bam";
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    print $pbs "
if [ ! -s $unsorted_bam_file ]; then
  echo bwa_mem=`date`
  bwa mem $option -R '$rg' $bwa_index $sample_files_str | samtools view -bS -o $unsorted_bam_file
  bwa 2>\&1 | grep Version | cut -d ' ' -f2 | cut -d '-' -f1 | awk '{print \"bwa,v\"\$1}' > ${final_file}.bwa.version
fi
";
    my $rmlist = "";
    if ($cleansam) {
      print $pbs "
if [[ -s $unsorted_bam_file && ! -s $clean_bam_file ]]; then
  echo CleanSam=`date`
  java -jar $picard_jar CleanSam VALIDATION_STRINGENCY=SILENT I=$unsorted_bam_file O=$clean_bam_file
fi
";
      $rmlist            = $rmlist . " " . $unsorted_bam_file;
      $unsorted_bam_file = $clean_bam_file;
    }

    if ($mark_duplicates) {
      print $pbs "
if [ -s $unsorted_bam_file && ! -s $rmdup_bam_file ]; then
  echo MarkDuplicate=`date` 
  java -jar $picard_jar MarkDuplicates \\
    INPUT=$unsorted_bam_file \\
    OUTPUT=$rmdup_bam_file \\
    METRICS_FILE=${final_file}.metrics \\
    VALIDATION_STRINGENCY=SILENT \\
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \\
    ASSUME_SORT_ORDER=\"queryname\" \\
    CLEAR_DT=\"false\" \\
    ADD_PG_TAG_TO_READS=false     
fi
";
      $rmlist            = $rmlist . " $bam_file ${bam_file}.bai";
      $unsorted_bam_file = $rmdup_bam_file;
    }

    if ($sortByCoordinate) {
      my $chromosome_grep_command = getChromosomeFilterCommand( $final_file, $chromosome_grep_pattern );

      print $pbs "    
if [ -s $unsorted_bam_file ]; then
  echo sort_bam=`date`
  samtools sort -@ $thread -m $sort_memory $unsorted_bam_file -o $final_file
  samtools index $final_file 
  $chromosome_grep_command
fi
";
      $rmlist = $rmlist . " " . $unsorted_bam_file;
    }
    else {
      print $pbs "
if [ -s $unsorted_bam_file ]; then
  mv $unsorted_bam_file $final_file
fi
";
    }

    print $pbs "
if [ -s $final_file ]; then
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

  print "!!!shell file $shfile created, you can run this shell file to submit all BWA tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $sortByCoordinate = get_option( $config, $section, "sort_by_coordinate", 1 );
  my $mark_duplicates  = hasMarkDuplicate( $config->{$section} );
  my %raw_files        = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $final_file = $sample_name . ( $mark_duplicates ? ".rmdup" : "" ) . ( $sortByCoordinate ? ".sorted" : ".unsorted" ) . ".bam";
    $final_file = "${result_dir}/$final_file";
    my @result_files = ();
    push( @result_files, $final_file );
    push( @result_files, $final_file . ".stat" );
    push( @result_files, $final_file . ".bwa.version" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
