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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $bowtie1_index = $config->{$section}{bowtie1_index} or die "define ${section}::bowtie1_index first";
  my $samformat  = get_option( $config, $section, "samformat",  1 );
  my $mappedonly = get_option( $config, $section, "mappedonly", 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . " \n";

  if ($samformat) {
    my $samonly = get_option( $config, $section, "samonly", 0 );
    my $sortbam = get_option( $config, $section, "sortbam", 1 );

    for my $sample_name ( sort keys %raw_files ) {
      my @sample_files = @{ $raw_files{$sample_name} };
      my $sam_file     = $sample_name . ".sam";

      my $bowtiesam        = $sam_file;
      my $mappedonlycmd    = "";
      my $mappedonlyoption = "";
      if ($mappedonly) {
        $bowtiesam     = $sample_name . ".all.sam";
        $mappedonlycmd = "
if [ -s $bowtiesam ]; then
  samtools view -F 4 $bowtiesam > $sam_file
  rm $bowtiesam
fi";
        $mappedonlyoption = "-F 4";
      }

      my $bam_file = $sample_name . ".bam";

      my $indent = "";
      my $tag    = "--sam-RG ID:$sample_name --sam-RG LB:$sample_name --sam-RG SM:$sample_name --sam-RG PL:ILLUMINA";

      my $bowtie1_aln_command;
      if ( $sample_files[0] =~ /.gz$/ ) {
        if ( scalar(@sample_files) == 1 ) {
          $bowtie1_aln_command = "zcat $sample_files[0] | bowtie $option -S $tag $bowtie1_index - $bowtiesam";
        }
        else {
          my $f1 = $sample_files[0];
          my $f2 = $sample_files[1];

          $bowtie1_aln_command = "
mkfifo ${f1}.fifo
zcat $f1 > ${f1}.fifo &

mkfifo ${f2}.fifo
zcat $f2 > ${f2}.fifo &
        
bowtie $option -S $tag $bowtie1_index ${f1}.fifo,${f2}.fifo $bowtiesam
 
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
      my $cur_dir  = create_directory_or_die( $result_dir . "/$sample_name" );

      print $sh "\$MYCMD ./$pbs_name \n";

      my $log_desc = $cluster->get_log_description($log);

      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir );

      if ($samonly) {
        print $pbs "
if [ -s $sam_file ]; then
  echo job has already been done. if you want to do again, delete $sam_file and submit job again.
  exit 0
fi

$bowtie1_aln_command

$mappedonlycmd
";
      }
      else {
        print $pbs "
if [ -s $bam_file ]; then
  echo job has already been done. if you want to do again, delete $bam_file and submit job again.
  exit 0
fi

$bowtie1_aln_command

if [ -s $bowtiesam ]; then
";
        if ($sortbam) {
          print $pbs "  samtools view -Shu $mappedonlyoption $bowtiesam | samtools sort -o $bam_file -
  if [ -s $bam_file ]; then
    samtools index $bam_file 
    samtools flagstat $bam_file > ${bam_file}.stat
";
        }
        else {
          print $pbs "samtools view -S $mappedonlyoption -b $bowtiesam > ${sample_name}.bam
  if [ -s $bam_file ]; then
";
        }
        print $pbs "    rm $bowtiesam
  fi
fi
";
      }
      
      $self->close_pbs($pbs);

      print "$pbs_file created\n";
    }
    close $sh;

    if ( is_linux() ) {
      chmod 0755, $shfile;
    }
  }
  else {
    for my $sample_name ( sort keys %raw_files ) {
      my @sample_files = @{ $raw_files{$sample_name} };
      my $final_file   = $sample_name . ".out";

      my $indent = "";

      my $bowtie1_aln_command;
      if ( $sample_files[0] =~ /.gz$/ ) {
        if ( scalar(@sample_files) == 1 ) {
          $bowtie1_aln_command = "zcat $sample_files[0] | bowtie $option $bowtie1_index - $final_file";
        }
        else {
          my $f1 = $sample_files[0];
          my $f2 = $sample_files[1];

          $bowtie1_aln_command = "
mkfifo ${f1}.fifo
zcat $f1 > ${f1}.fifo &

mkfifo ${f2}.fifo
zcat $f2 > ${f2}.fifo &
        
bowtie $option $bowtie1_index ${f1}.fifo,${f2}.fifo $final_file
 
rm ${f1}.fifo
rm ${f2}.fifo
";
        }
      }
      else {
        my $fastqs = join( ',', @sample_files );
        $bowtie1_aln_command = "bowtie $option $bowtie1_index $fastqs $final_file ";
      }

      my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
      my $pbs_name = basename($pbs_file);
      my $log      = $self->get_log_filename( $log_dir, $sample_name );
      my $log_desc = $cluster->get_log_description($log);

      my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

      print $sh "\$MYCMD ./$pbs_name \n";

      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

      print $pbs "
$bowtie1_aln_command
";
      $self->close_pbs($pbs);
      print "$pbs_file created\n";
    }
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

1;
