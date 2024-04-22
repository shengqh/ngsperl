#!/usr/bin/perl
package GATK4::ApplyBQSRScatter;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use GATK4::IntervalsScatterTask;

our @ISA = qw(GATK4::IntervalsScatterTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_ab";
  $self->{_docker_prefix}   = "gatk4_";
  bless $self, $class;

  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $scatter_map = get_interval_file_map($config, $section);

  my $expect_bam_suffix = get_option( $config, $section, "bam_suffix", ".bam" );
  my $expect_bam_index_suffix = get_option( $config, $section, "bam_index_suffix", ".bai" );

  my $ref_fasta = get_option_file( $config, $section, "fasta_file", 1 );

  my $java_option = $self->get_java_option($config, $section, $memory);

  #$self->get_docker_value(0);

  my %bam_files = %{ get_raw_files( $config, $section ) };
  my $bqsr_report_files = get_raw_files( $config, $section, "bqsr_report_files" );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %bam_files ) {
    my @sample_files = @{ $bam_files{$sample_name} };
    my $input_bam     = $sample_files[0];

    my ($bam_suffix) = $input_bam =~ /(\.[^.]+)$/;
    my $bam_index_suffix = $bam_suffix eq ".cram" ? ".crai" : ".bai";

    die("expecting suffix $expect_bam_suffix, but found $bam_suffix") if $expect_bam_suffix ne $bam_suffix;

    my $recalibration_report = $bqsr_report_files->{$sample_name}[0];

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    for my $scatter_name (sort keys %$scatter_map) {
      my $interval_file = $scatter_map->{$scatter_name};
      my $prefix = get_key_name($sample_name, $scatter_name);
      my $final_file = $prefix . ".recalibrated$bam_suffix";
      my $final_index = $prefix . ".recalibrated$bam_index_suffix";
      
      my $pbs_file = $self->get_pbs_filename( $pbs_dir, $prefix );
      my $pbs_name = basename($pbs_file);
      my $log      = $self->get_log_filename( $log_dir, $prefix );

      print $sh "if [[ ! -s $cur_dir/$final_index ]]; then 
  \$MYCMD ./$pbs_name 
fi
";

      my $log_desc = $cluster->get_log_description($log);
      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_index );
      print $pbs "
gatk --java-options \"$java_option\" \\
  ApplyBQSR \\
  --add-output-sam-program-record \\
  -R ${ref_fasta} \\
  -I ${input_bam} \\
  --use-original-qualities \\
  -O $final_file \\
  -bqsr ${recalibration_report} \\
  -L ${interval_file}

status=\$?
if [[ \$status -ne 0 ]]; then
  touch $prefix.failed
  rm -f $prefix.succeed $final_file $final_index
else
  touch $prefix.succeed
  rm -f $prefix.failed
fi
";
      
      $self->close_pbs( $pbs, $pbs_file );
    }
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks \n ";
}

sub get_result_files {
  my ( $self, $config, $section, $result_dir, $sample_name, $scatter_name, $key_name ) = @_;

  my $expect_bam_suffix = get_option( $config, $section, "bam_suffix", ".bam" );

  my $final_file = "${result_dir}/${sample_name}/${key_name}.recalibrated$expect_bam_suffix";
  return [$final_file];
}

1;
