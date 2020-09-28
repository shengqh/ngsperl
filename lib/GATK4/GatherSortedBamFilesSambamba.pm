#!/usr/bin/perl
package GATK4::GatherSortedBamFilesSambamba;

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
  $self->{_docker_prefix}   = "cqs_";
  bless $self, $class;

  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $java_option = $self->get_java_option($config, $section, $memory);

  $self->get_docker_value(0);

  my $extension = get_option($config, $section, "extension", ".bam");

  my $gather_map = $self->get_gather_map($config, $section);
  my $raw_files = get_raw_files( $config, $section );
  my $gather_file_map = $self->get_gather_file_map($gather_map, $raw_files);

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $gather_name (sort keys %$gather_file_map){
    my $key_files = $gather_file_map->{$gather_name};
    my $input_bam_files = join(" \\\n  ", @$key_files);

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $gather_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $gather_name );

    my $final_file = $gather_name . $extension;
    my $final_index = $final_file . ".bai";

    print $sh "if [[ ! -s $result_dir/$final_file ]]; then 
  \$MYCMD ./$pbs_name 
fi
";

    my $tmp_file = $gather_name . ".tmp" . $extension;
    my $tmp_index = $tmp_file . ".bai";

    my $log_desc = $cluster->get_log_description($log);
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_index );
    print $pbs "
sambamba merge $option -t $thread $tmp_file \\
  $input_bam_files

if [[ -s $tmp_index ]]; then
  mv $tmp_file $final_file
  mv $tmp_index $final_index
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
  my $extension = get_option($config, $section, "extension", ".bam");
  my $final_file = "${result_dir}/${gather_name}${extension}";
  return [$final_file];
}

1;
