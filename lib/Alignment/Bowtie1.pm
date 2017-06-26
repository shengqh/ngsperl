#!/usr/bin/perl
package Alignment::Bowtie1;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use Alignment::AbstractBowtie;
use Alignment::AlignmentUtils;

our @ISA = qw(Alignment::AbstractBowtie);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_bt1";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $bowtie1_index = $config->{$section}{bowtie1_index} or die "define ${section}::bowtie1_index first";
  my $sortbam                 = get_option( $config, $section, "sortbam",                 1 );
  my $mappedonly              = get_option( $config, $section, "mappedonly",              0 );
  my $chromosome_grep_pattern = get_option( $config, $section, "chromosome_grep_pattern", "" );
  my $outputToSameFolder      = get_option( $config, $section, "output_to_same_folder",   1 );

  my $mark_duplicates         = hasMarkDuplicate( $config->{$section} );
  my $picard_jar = "";
  if($mark_duplicates){
    $picard_jar = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1 );
  }

  $option = $option . " -p $thread";

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . " \n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sam_file     = $sample_name . ".sam";

    my $bowtiesam        = $sam_file;
    my $alignlog         = $sample_name . ".log";
    my $mappedonlyoption = "";
    if ($mappedonly) {
      $bowtiesam        = $sample_name . ".all.sam";
      $mappedonlyoption = "-F 4";
    }

    my $indent = "";
    my $tag    = "--sam-RG ID:$sample_name --sam-RG LB:$sample_name --sam-RG SM:$sample_name --sam-RG PL:ILLUMINA --sam-RG PU:$sample_name";

    my $bowtie1_aln_command;
    if ( $sample_files[0] =~ /.gz$/ ) {
      if ( scalar(@sample_files) == 1 ) {
        $bowtie1_aln_command = "zcat $sample_files[0] | bowtie $option -S $tag $bowtie1_index - $bowtiesam 2>$alignlog";
      }
      else {
        my $f1 = $sample_files[0];
        my $f2 = $sample_files[1];

        $bowtie1_aln_command = "mkfifo ${f1}.fifo
  zcat $f1 > ${f1}.fifo &

  mkfifo ${f2}.fifo
  zcat $f2 > ${f2}.fifo &
        
  bowtie $option -S $tag $bowtie1_index ${f1}.fifo,${f2}.fifo $bowtiesam 2>$alignlog
 
  rm ${f1}.fifo
  rm ${f2}.fifo";
      }
    }
    else {
      my $fastqs = join( ',', @sample_files );
      $bowtie1_aln_command = "bowtie $option -S $tag $bowtie1_index $fastqs $bowtiesam";
    }

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $cur_dir  = $outputToSameFolder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );

    my $bam_file = $sample_name . ".bam";
    my $final_file = $mark_duplicates ? $sample_name . ".rmdup.bam" : $bam_file;

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    print $pbs "
if [[ ! -s $bam_file && ! -s $bowtiesam ]]; then
  $bowtie1_aln_command 
fi
";
    if ($sortbam) {
      my $chromosome_grep_command = getChromosomeFilterCommand( $bam_file, $chromosome_grep_pattern );

      print $pbs "
if [[ -s $bowtiesam && ! -s $bam_file ]]; then
  samtools view -Shu $mappedonlyoption $bowtiesam | samtools sort -@ $thread -T ${sample_name}_tmp -o $bam_file -
  if [ -s $bam_file ]; then
    rm $bowtiesam
    samtools index $bam_file 
    $chromosome_grep_command
  fi
fi
";
      if ($mark_duplicates) {
        print $pbs "
if [ -s $bam_file ]; then
  echo MarkDuplicate=`date` 
  java -jar $picard_jar MarkDuplicates I=$bam_file O=$final_file ASSUME_SORTED=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT M=${final_file}.metrics
  if [ -s $final_file ]; then
    rm $bam_file ${bam_file}.bai
    samtools index $final_file 
  fi
fi
";
      }
    }
    else {
      print $pbs "
if [ -s $bowtiesam ]; then
  samtools view -S $mappedonlyoption -b $bowtiesam > ${sample_name}.bam
  if [ -s $bam_file ]; then
    rm $bowtiesam
  fi
fi
";
    }

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

1;
