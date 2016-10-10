#!/usr/bin/perl
package CQS::CQSDatatable;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_tb";
  bless $self, $class;
  return $self;
}

sub get_result {
  my ( $task_name, $option ) = @_;

  my $result;
  if ( $option =~ /-o\s+(\S+)/ ) {
    $result = $1;
  }
  else {
    $result = $task_name . ".count";
    $option = $option . " -o " . $result;
  }
  return ( $result, $option );
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $cqstools = get_cqstools( $config, $section, 1 );
  my $mapFile = get_param_file( $config->{$section}{name_map_file}, "name_map_file", 0 );

  my $mapoption = "";
  if ( defined $mapFile ) {
    $mapoption = "-m $mapFile";
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $filelist = $self->get_file( $pbs_dir, $task_name, ".filelist" );
  open( FL, ">$filelist" ) or die "Cannot create $filelist";
  for my $sample_name ( sort keys %raw_files ) {
    my @bam_files = @{ $raw_files{$sample_name} };
    my $bam_file  = $bam_files[0];
    print FL $sample_name, "\t", $bam_file, "\n";
  }
  close(FL);

  my ( $result_file, $newoption ) = get_result( $task_name, $option );

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $result_file );

  print $pbs "mono $cqstools data_table $newoption -l $filelist $mapoption";

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $result = {};
  my ( $result_file, $newoption ) = get_result( $task_name, $option );
  $result_file = $result_dir . "/" . $result_file;

  my $filelist = $pbs_dir . "/" . $self->get_name( $task_name, ".filelist" );

  my @result_files = ();
  push( @result_files, $result_file );
  push( @result_files, $filelist );

  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
