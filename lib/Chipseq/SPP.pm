#!/usr/bin/perl
package Chipseq::SPP;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::GroupTask;
use CQS::NGSCommon;
use Data::Dumper;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_spp";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $spp_r = get_directory( $config, $section, "spp_r", 1 );
  if ( $option eq "" ) {
    $option = "-npeak=300000 -savr -savp -savd -rf ";
  }

  my %treatments_files = %{ $self->get_grouped_raw_files( $config, $section, "groups" ) };

  my %control_files;
  if ( has_raw_files( $config, $section, "inputs" ) ) {
    %control_files = %{ $self->get_grouped_raw_files( $config, $section, "inputs" ) };
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $group_name ( sort keys %treatments_files ) {
    my @sample_files = @{ $treatments_files{$group_name} };
    my $treatment    = "-c=" . $sample_files[0];

    my $control = "";
    if (%control_files) {
      my @c_files = @{ $control_files{$group_name} };
      $control = "-i=" . $c_files[0];
    }

    my $cur_dir = create_directory_or_die( $result_dir . "/$group_name" );
    my $tmp_dir = create_directory_or_die( $result_dir . "/${group_name}/tmp" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $group_name );

    my $log_desc   = $cluster->get_log_description($log);
    my $final_file = "${cur_dir}/${group_name}_spp.tab";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );
    print $pbs "Rscript $spp_r $option $treatment $control -odir=$cur_dir -out=${group_name}_spp.tab -tmpdir=$tmp_dir \n";
    $self->close_pbs( $pbs, $pbs_file );

    print $sh "\$MYCMD ./$pbs_name \n";
  }
  print $sh "exit 0\n";
  close $sh;

  print "!!!shell file $shfile created, you can run this shell file to submit tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %treatments_files = %{ $self->get_grouped_raw_files( $config, $section, "groups" ) };

  my $result = {};
  for my $group_name ( sort keys %treatments_files ) {
    my $cur_dir      = $result_dir . "/$group_name";
    my @result_files = ();
    my $final_file   = "${cur_dir}/${group_name}_spp.tab";
    push( @result_files, $final_file );
    $result->{$group_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
