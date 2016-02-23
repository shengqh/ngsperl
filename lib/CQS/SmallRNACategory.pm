#!/usr/bin/perl
package CQS::SmallRNACategory;

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
  $self->{_suffix} = "_cat";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $cqstools = get_cqstools( $config, $section, 1 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $pdfoption = get_option( $config, $section, "pdfgraph", 0 ) ? "--pdf" : "";

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  my $filelist = $self->get_file( $pbs_dir, $task_name, ".filelist" );
  open( my $fl, ">$filelist" ) or die "Cannot create $filelist";

  if ( defined $config->{$section}{groups} || defined $config->{$section}{groups_ref} ) {
    my $groups = get_raw_files( $config, $section, "groups" );
    for my $group_name ( sort keys %{$groups} ) {
      my @samples = @{ $groups->{$group_name} };
      for my $sample_name ( sort @samples ) {
        my @smallRNAFiles = @{ $raw_files{$sample_name} };
        my $smallRNAFile  = $smallRNAFiles[0];

        print $fl $group_name, "\t", $sample_name, "\t", $smallRNAFile, "\n";
      }
    }
  }
  else {
    for my $sample_name ( sort keys %raw_files ) {
      my @smallRNAFiles = @{ $raw_files{$sample_name} };
      my $smallRNAFile  = $smallRNAFiles[0];

      print $fl $task_name, "\t", $sample_name, "\t", $smallRNAFile, "\n";
    }
  }
  close($fl);

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, "${task_name}.catcount" );
  print $pbs "mono $cqstools smallrna_group -i $filelist -o $result_dir";
  $self->close_pbs( $pbs, $pbs_file );

  print $sh "\$MYCMD ./$pbs_name \n";
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all smallRNA category tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my @result_files = ("${result_dir}/${task_name}.catcount");
  my $result       = {};
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
