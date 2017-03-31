#!/usr/bin/perl
package Alignment::Bowtie2;

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
  $self->{_suffix} = "_bt2";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $bowtie2_index = $config->{$section}{bowtie2_index} or die "define ${section}::bowtie2_index first";
  my $chromosome_grep_pattern = get_option( $config, $section, "chromosome_grep_pattern", "" );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct), "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sam_file     = $sample_name . ".sam";
    my $bam_file     = $sample_name . ".bam";

    my $indent = "";
    my $tag    = "--sam-RG ID:$sample_name --sam-RG LB:$sample_name --sam-RG SM:$sample_name --sam-RG PL:ILLUMINA";

    my $input = "";
    if ( scalar(@sample_files) == 2 ) {    #pairend data
      $input = "-1 " . $sample_files[0] . " -2 " . $sample_files[1];
    }
    else {
      my $fastqs = join( ',', @sample_files );
      $input = "-u " . $fastqs;
    }
    my $bowtie2_aln_command = "bowtie2 -p $thread $option -x $bowtie2_index $input $tag -S $sam_file";

    my $index_command = get_index_command( $bam_file, $indent );
    my $stat_command = get_stat_command( $bam_file, $indent );

    my $pbs_name = $self->pbs_name($sample_name);
    my $pbs_file = $pbs_dir . "/$pbs_name";
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    my $final_file              = $bam_file;
    my $chromosome_grep_command = getChromosomeFilterCommand( $bam_file, $chromosome_grep_pattern );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $bam_file );

    print $pbs "

$bowtie2_aln_command

if [ -s $sam_file ]; then
  samtools view -Shu -F 256 $sam_file | samtools sort -o $bam_file -
  if [ -s $bam_file ]; then
    samtools index $bam_file 
    $chromosome_grep_command
    samtools flagstat $bam_file > ${bam_file}.stat 
    rm $sam_file
  fi
fi
";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

1;
