#!/usr/bin/perl
package CQS::Bam2Fastq;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_b2q";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $ispaired = get_option( $config, $section, "ispaired", 0 );

  my $sort_before_convert = get_option( $config, $section, "sort_before_convert" );
  my $sort_thread         = get_option( $config, $section, "sort_thread" );
  my $sortoption = "";

  my $unmapped_only = get_option( $config, $section, "unmapped_only", 0 );
  my $unzipped      = get_option( $config, $section, "unzipped",      0 );
  if ($unzipped) {
    $option = $option . " -u";
  }

  print "sort_before_convert = $sort_before_convert\n";
  print "sort_thread = $sort_thread\n";

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $bam_file     = $sample_files[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $final_file = $ispaired ? $sample_name . ".1.fastq" : $sample_name . ".fastq";
    if ( !$unzipped ) {
      $final_file = $final_file . ".gz";
    }
    my $convertCmd = "";

    if ($sort_before_convert) {
      my $sourceFile = "${sample_name}.sortname.bam";
      my $sortCmd;
      if ($unmapped_only) {
        $sortCmd = "samtools view -b -f 4 $bam_file | samtools sort $option -n $sortoption - ${sample_name}.sortname";
      }
      else {
        $sortCmd = "samtools sort $option -n $sortoption $bam_file ${sample_name}.sortname";
      }

      $convertCmd = "if [ ! -s $sourceFile ]; then
    $sortCmd
  fi
  
  cqstools bam2fastq $option -i $sourceFile -o $sample_name 
  
  if [ -s $final_file ]; then
    rm $sourceFile
  fi";
    }
    else {
      if ($unmapped_only) {
        my $unmapped_bam = "${sample_name}.unmapped.bam";
        $convertCmd = "samtools view -b -f 4 $bam_file > $unmapped_bam
  cqstools bam2fastq $option -i $unmapped_bam -o $sample_name
  rm $unmapped_bam ";
      }
      else {
        $convertCmd = "cqstools bam2fastq $option -i $bam_file -o $sample_name ";
      }
    }

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    print $pbs $convertCmd;

    $self->close_pbs( $pbs, $pbs_file );

  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all Bam2Fastq tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $ispaired = get_option( $config, $section, "ispaired", 0 );
  my $unzipped = get_option( $config, $section, "unzipped", 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {

    my @result_files = ();
    if ($ispaired) {
      if ($unzipped) {
        push( @result_files, $result_dir . "/" . $sample_name . ".1.fastq" );
        push( @result_files, $result_dir . "/" . $sample_name . ".2.fastq" );
      }
      else {
        push( @result_files, $result_dir . "/" . $sample_name . ".1.fastq.gz" );
        push( @result_files, $result_dir . "/" . $sample_name . ".2.fastq.gz" );
      }
    }
    else {
      if ($unzipped) {
        push( @result_files, $result_dir . "/" . $sample_name . ".fastq" );
      }
      else {
        push( @result_files, $result_dir . "/" . $sample_name . ".fastq.gz" );
      }
    }

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
