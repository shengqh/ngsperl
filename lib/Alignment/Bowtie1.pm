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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $bowtie1_index = $config->{$section}{bowtie1_index} or die "define ${section}::bowtie1_index first";
  my $samformat               = get_option( $config, $section, "samformat",               1 );
  my $mappedonly              = get_option( $config, $section, "mappedonly",              0 );
  my $chromosome_grep_pattern = get_option( $config, $section, "chromosome_grep_pattern", "" );
  my $outputToSameFolder      = get_option( $config, $section, "output_to_same_folder",   0 );
  my $add_RG_to_read = get_option( $config, $section, "add_RG_to_read", 0 );
  my $picard;
  if($add_RG_to_read){
    $picard = get_param_file($config->{$section}{"picard_jar"}, "picard_jar", 1);
  }

  $option = $option . " -p $thread";

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . " \n";

  if ($samformat) {
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
      my $tag    = "--sam-RG ID:$sample_name --sam-RG LB:$sample_name --sam-RG SM:$sample_name --sam-RG PL:ILLUMINA --sam-RG PU:$sample_name";

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
      my $cur_dir  = $outputToSameFolder ? $result_dir : create_directory_or_die( $result_dir . "/$sample_name" );

      my $chromosome_grep_command = "";
      my $final_file              = $bam_file;
      if ( $sortbam && ( $chromosome_grep_pattern ne "" ) ) {
        my $tmp_file = $sample_name . ".filtered.bam";
        $chromosome_grep_command = "
    echo filtering bam by chromosome pattern $chromosome_grep_pattern
    samtools idxstats $bam_file | cut -f 1 | grep $chromosome_grep_pattern | xargs samtools view -b $bam_file > $tmp_file
    samtools flagstat $bam_file > ${bam_file}.raw.stat
    rm $bam_file
    rm ${bam_file}.bai
    mv $tmp_file $bam_file
    samtools index $bam_file
";
      }
      
      my $add_RG_to_read_command = "";
      if($add_RG_to_read){
        my $tmp_file = $sample_name . ".rg.bam";
        $add_RG_to_read_command = "
    echo add_RG_to_read_command by picard
    java -jar $picard AddOrReplaceReadGroups I=$bam_file O=$tmp_file ID=$sample_name LB=$sample_name SM=$sample_name PL=illumina PU=$sample_name CREATE_INDEX=False
    rm $bam_file
    rm ${bam_file}.bai
    mv $tmp_file $bam_file
    samtools index $bam_file
";
      }

      print $sh "\$MYCMD ./$pbs_name \n";

      my $log_desc = $cluster->get_log_description($log);

      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

      print $pbs "$bowtie1_aln_command \n";
      if ($sortbam) {
        print $pbs "
if [ -s $bowtiesam ]; then
  samtools view -Shu $mappedonlyoption $bowtiesam | samtools sort -@ $thread -T ${sample_name}_tmp -o $bam_file -
  if [ -s $bam_file ]; then
    rm $bowtiesam
    samtools index $bam_file 
    $chromosome_grep_command
    $add_RG_to_read_command    
    samtools flagstat $bam_file > ${bam_file}.stat
  fi
fi
";
      }
      else {
        print $pbs "
if [ -s $bowtiesam ]; then
  samtools view -S $mappedonlyoption -b $bowtiesam > ${sample_name}.bam
  if [ -s $bam_file ]; then
    rm $bowtiesam
    $add_RG_to_read_command    
    samtools flagstat $bam_file > ${bam_file}.stat
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
