#!/usr/bin/perl
package GATK4::GenomicsDBImportScatter;

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
  $self->{_suffix} = "_gbi";
  $self->{_docker_prefix}   = "gatk4_";
  bless $self, $class;

  return $self;
}

sub get_sample_names {
  my ($self, $config, $section) = @_;
  my ( $task_name ) = $self->init_parameter( $config, $section, 0 );
  return [$task_name];
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $scatter_map = get_interval_file_map($config, $section);

  my $java_option = $self->get_java_option($config, $section, $memory);

  my $gvcf_files = get_raw_files( $config, $section );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  my $sample_name_map_file = $self->get_file( $result_dir, $task_name, ".sample_name_map.txt" );
  open( my $sf, ">$sample_name_map_file" ) or die "Cannot create $sample_name_map_file";
  for my $sample_name ( sort keys %$gvcf_files ) {
    my $gvcf_file = $gvcf_files->{$sample_name}[0];
    print $sf $sample_name . "\t" . $gvcf_file . "\n";
  }
  close($sf);

  for my $scatter_name (sort keys %$scatter_map) {
    my $interval_file = $scatter_map->{$scatter_name};
    my $prefix = get_key_name($task_name, $scatter_name);
    
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $prefix );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $prefix );

    my $final_folder = "$result_dir/$prefix";

    print $sh "if [[ ! -e $final_folder ]]; then 
\$MYCMD ./$pbs_name 
fi
";

    my $log_desc = $cluster->get_log_description($log);
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_folder );
    print $pbs "
rm -rf $prefix

gatk --java-options \"-Xmx8g -Xms8g\"  \\
  GenomicsDBImport \\
  --genomicsdb-workspace-path $prefix \\
  --batch-size 50 \\
  -L $interval_file \\
  --sample-name-map $sample_name_map_file \\
  --reader-threads 5 \\
  --merge-input-intervals \\
  --consolidate

";
    
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks \n ";
}

sub get_result_files {
  my ( $self, $config, $section, $result_dir, $sample_name, $scatter_name, $key_name ) = @_;
  my $final_folder = "${result_dir}/${key_name}";
  return [$final_folder];
}

1;
