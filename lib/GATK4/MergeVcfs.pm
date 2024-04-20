#!/usr/bin/perl
package GATK4::MergeVcfs;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use GATK4::GATK4UniqueTask;

our @ISA = qw(GATK4::GATK4UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_mv";
  $self->{_docker_prefix}   = "gotc_";
  bless $self, $class;

  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $java_option = $self->get_java_option($config, $section, $memory);

  $self->get_docker_value(1);

  my $raw_files = get_raw_files( $config, $section );

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $filelist = $self->get_file($result_dir, $task_name, ".input.list");
  writeFileList($filelist, $raw_files, 1, 1, 1);

  my $filelist_basename = basename($filelist);

  my $extension = get_option($config, $section, "extension", 1);
  my $final_file = $task_name . $extension;
  my $final_index = $final_file . ".tbi";

  my $tmp_file = $task_name . ".tmp" . $extension;
  my $tmp_index = $tmp_file . ".tbi";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_index );
  print $pbs "
gatk --java-options \"-Xms${memory} -Xmx${memory}\" \\
      MergeVcfs \\
      -I $filelist_basename \\
      -O $tmp_file

status=\$?
if [[ \$status -eq 0 ]]; then
  touch ${final_file}.succeed
  rm -f ${final_file}.failed
  mv $tmp_file $final_file
  mv $tmp_index $final_index
else
  touch ${final_file}.failed
  rm -f ${final_file}.succeed
  rm -f $tmp_file
  rm -f $tmp_index
fi
";
    
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  my $result = {};

  my $extension = get_option($config, $section, "extension", 1);
  my $final_file = $task_name . $extension;
  $result->{$task_name} = filter_array( ["$result_dir/$final_file"], $pattern );
  return $result;
}

1;
