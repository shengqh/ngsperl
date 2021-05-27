#!/usr/bin/perl
package CQS::NGSCommon;

use strict;
use warnings;

use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(should_use_tmp_folder get_bam_tag get_sorted_bam get_sam2bam_command get_sort_index_command get_sort_command get_index_command get_stat_command transcript_gtf_index_exists)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

#if target folder is at local storage, tmp folder should not be used
sub should_use_tmp_folder {
  my ($target_dir) = @_;
  if($target_dir =~ /^\/workspace/ || $target_dir =~ /^\/data1/){
    return(0)
  }else{
    return(1);
  }
}

sub get_bam_tag {
  my ($sample_name) = @_;
  my $tag = "'\@RG ID:$sample_name LB:$sample_name SM:$sample_name PL:ILLUMINA'";
  return ($tag);
}

sub get_sorted_bam {
  my $bam_file         = shift;
  my $bamSortedPrefix = change_extension( $bam_file, "_sorted" );
  my $bamSortedFile   = $bamSortedPrefix . ".bam";
  return ( $bamSortedFile, $bamSortedPrefix );
}

sub get_sam2bam_command {
  my ( $sam_file, $bam_file, $indent ) = @_;

  if ( !defined($indent) ) {
    $indent = "";
  }

  my $command = "${indent}if [ -s $sam_file ]; then
${indent}  echo sam2bam=`date`
${indent}  samtools view -b -S $sam_file -o $bam_file
${indent}fi";

  return ($command);
}

sub get_sort_command {
  my ( $bam_file, $bamSortedPrefix, $indent ) = @_;

  if ( !defined($indent) ) {
    $indent = "";
  }

  my $bamSortedFile;
  if ( !defined $bamSortedPrefix ) {
    ( $bamSortedFile, $bamSortedPrefix ) = get_sorted_bam($bam_file);
  }
  else {
    $bamSortedFile = $bamSortedPrefix . ".bam";
  }

  my $command = "${indent}if [ -s $bam_file ]; then
${indent}  echo BamSort=`date` 
${indent}  samtools sort $bam_file $bamSortedPrefix 
${indent}fi";
  return ($command);
}

sub get_index_command {
  my ( $bamSortedFile, $indent ) = @_;

  if ( !defined($indent) ) {
    $indent = "";
  }

  my $command = "${indent}if [ -s $bamSortedFile ]; then
${indent}  echo BamIndex=`date` 
${indent}  samtools index $bamSortedFile
${indent}fi";
  return ($command);
}

sub get_stat_command {
  my ( $bamSortedFile, $indent ) = @_;

  if ( !defined($indent) ) {
    $indent = "";
  }

  my $command = "${indent}if [ -s $bamSortedFile ]; then
${indent}  echo BamStat=`date`
${indent}  samtools flagstat $bamSortedFile > ${bamSortedFile}.stat 
${indent}fi";

  return ($command);
}

sub get_sort_index_command {
  my ( $bam_file, $bamSortedPrefix, $indent ) = @_;
  my $bamSortedFile;
  if ( !defined $bamSortedPrefix ) {
    ( $bamSortedFile, $bamSortedPrefix ) = get_sorted_bam($bam_file);
  }
  else {
    $bamSortedFile = $bamSortedPrefix . ".bam";
  }

  return get_sort_command($bam_file, $bamSortedPrefix, $indent) . "\n" . get_index_command($bamSortedFile, $indent);  
}

sub transcript_gtf_index_exists {
  my $transcript_gtf_index = shift;
  my $result               = 0;
  if ( defined($transcript_gtf_index) ) {
    my $file = $transcript_gtf_index . ".rev.1.bt2";
    $result = -e $file;
  }
  return ($result);
}

1;
