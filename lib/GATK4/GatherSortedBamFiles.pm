#!/usr/bin/perl
package GATK4::GatherSortedBamFiles;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use GATK4::IntervalsGatherTask;

our @ISA = qw(GATK4::IntervalsGatherTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_gbf";
  $self->{_docker_prefix}   = "gatk4_";
  bless $self, $class;

  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $java_option = $self->get_java_option($config, $section, $memory);

  $self->get_docker_value(1);

  my $ref_fasta = get_option_file( $config, $section, "ref_fasta");

  my $gather_map = $self->get_gather_map($config, $section);
  my $raw_files = get_raw_files( $config, $section );
  my $gather_file_map = $self->get_gather_file_map($gather_map, $raw_files);

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $gather_name (sort keys %$gather_file_map){
    my $key_files = $gather_file_map->{$gather_name};
    my $input_bam_files = join(" \\\n  --INPUT ", @$key_files);

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $gather_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $gather_name );

    my $final_file = $gather_name . ".bam";
    my $final_index = $gather_name . ".bai";

    print $sh "
if [[ ! -s $final_file ]]; then
  \$MYCMD ./$pbs_name
fi
";

    my $log_desc = $cluster->get_log_description($log);
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
    print $pbs "
gatk --java-options \"$java_option\" \\
  GatherBamFiles \\
  --INPUT ${input_bam_files} \\
  --OUTPUT ${final_file} \\
  --CREATE_INDEX true \\
  --CREATE_MD5_FILE false

status=\$?
if [[ \$status -ne 0 ]]; then
  touch $final_file.failed
  rm -f $final_file.succeed $final_file $final_index 
else
  touch $final_file.succeed
  rm -f $final_file.failed
fi
";
    
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks. \n ";
}

sub get_result_files {
  my ( $self, $config, $section, $result_dir, $gather_name ) = @_;
  my $final_file = "${result_dir}/${gather_name}.bam";
  return [$final_file];
}

1;
