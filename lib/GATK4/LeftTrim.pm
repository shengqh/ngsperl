#!/usr/bin/perl
package GATK4::LeftTrim;

use strict;
use warnings;
use File::Basename;
use List::Util qw[min];
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use GATK4::GATK4Task;
use GATK4::VariantFilterUtils;

our @ISA = qw(GATK4::GATK4Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_lf";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  $self->get_docker_value(1);

  my $extension = get_option($config, $section, "extension", 1);

  my $vcf_files = get_raw_files( $config, $section );

  my $script = dirname(__FILE__) . "/fixLeftTrimDeletion.py";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $interval_name (sort keys %$vcf_files) {
    my $pass_file = $vcf_files->{$interval_name}[0];
    my $final_file = $interval_name . $extension;

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $interval_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $interval_name );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "if [[ ! -s $result_dir/$final_file ]]; then
  \$MYCMD ./$pbs_name 
fi
    
";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command );

    my ($left_trim_rmlist) = add_left_trim_pbs($self, $config, $section, $pbs, $interval_name, $pass_file, $final_file);

    print $pbs "
  if [[ -s $final_file ]]; then
    rm $left_trim_rmlist
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

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $extension = get_option($config, $section, "extension", 1);
  my $result = {};
  my $vcf_files = get_raw_files( $config, $section );
  for my $interval_name (sort keys %$vcf_files) {
    my $final_file = $interval_name . $extension;
    $result->{$interval_name} = filter_array( ["$result_dir/$final_file"], $pattern );
  }

  return $result;
}

1;
