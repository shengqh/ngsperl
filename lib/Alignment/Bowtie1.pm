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
  $self->{log_result} = 1;
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $bowtie1_index = $config->{$section}{bowtie1_index} or die "define ${section}::bowtie1_index first";
  my $mappedonly = get_option( $config, $section, "mappedonly", 0 );
  my $chromosome_grep_pattern = get_option( $config, $section, "chromosome_grep_pattern", "" );
  my $outputToSameFolder = $self->getOutputToSameFolder( $config, $section );
  my $output_sort_by_coordinate = getSortByCoordinate( $config, $section );

  my $mark_duplicates = hasMarkDuplicate( $config->{$section} );
  my $picard_jar      = "";
  if ($mark_duplicates) {
    $picard_jar = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1, not $self->using_docker() );
  }

  my $add_RG_to_read = get_option( $config, $section, "add_RG_to_read", 0 );
  if ($add_RG_to_read) {
    $picard_jar = get_param_file( $config->{$section}{picard_jar}, "picard_jar", 1, not $self->using_docker() );
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

    my $bam_file = $sample_name . ".bam";
    my $final_file = ( $output_sort_by_coordinate && $mark_duplicates ) ? $sample_name . ".rmdup.bam" : $bam_file;

    my $m_option = ($option =~ /\-m/)? "--max ${final_file}.max.txt":""; 

    my $indent = "";
    my $tag    = "--sam-RG ID:$sample_name --sam-RG LB:$sample_name --sam-RG SM:$sample_name --sam-RG PL:ILLUMINA --sam-RG PU:$sample_name";

    my $fastqs = join( ',', @sample_files );
    my $bowtie1_aln_command = "bowtie $option $m_option -S $tag -x $bowtie1_index $fastqs $bowtiesam 2>$alignlog";

    my $cmd_file_exists = check_file_exists_command(\@sample_files, "  ");

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $cur_dir  = $outputToSameFolder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    print $pbs "
if [[ ! -s $bam_file && ! -s $bowtiesam ]]; then
  $cmd_file_exists
  $bowtie1_aln_command 
  bowtie --version | grep bowtie | grep version | cut -d ' ' -f3 | awk '{print \"bowtie,v\"\$1}' > ${final_file}.version
fi
";

    my $addRgCommand = getAddRgCommand( $picard_jar, $add_RG_to_read, $bam_file, $sample_name );

    if ($output_sort_by_coordinate) {
      my $chromosome_grep_command = getChromosomeFilterCommand( $bam_file, $chromosome_grep_pattern );

      print $pbs "
if [[ -s $bowtiesam && ! -s $bam_file ]]; then
  samtools view -Shu $mappedonlyoption $bowtiesam | samtools sort -T ${sample_name}_tmp -o $bam_file -
  if [ -s $bam_file ]; then
    rm $bowtiesam
    samtools index $bam_file 
    $chromosome_grep_command
    $addRgCommand
    samtools idxstats $bam_file > ${bam_file}.chromosome.count
  fi
fi
";
      if ($mark_duplicates) {
        print $pbs "
if [ -s $bam_file ]; then
  echo MarkDuplicate=`date` 
  java -jar $picard_jar MarkDuplicates I=$bam_file O=$final_file ASSUME_SORTED=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT M=${final_file}.metrics
  if [ -s $final_file ]; then
    rm $bam_file ${bam_file}.bai ${bam_file}.chromosome.count
    samtools index $final_file 
    samtools idxstats $final_file > ${final_file}.chromosome.count
  fi
fi
";
      }
    }
    else {
      print $pbs "
if [ -s $bowtiesam ]; then
  samtools view -S $mappedonlyoption -b $bowtiesam > ${sample_name}.bam
  samtools idxstats ${sample_name}.bam > ${sample_name}.bam.chromosome.count
  if [ -s $bam_file ]; then
    $addRgCommand
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

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

1;
