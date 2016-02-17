#!/usr/bin/perl
package Blast::Blastn;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_bn";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $maximum_size = get_option( $config, $section, "maximum_size", 100 );
  my $bin_size     = get_option( $config, $section, "bin_size",     10 );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $blastn = dirname(__FILE__) . "/blastn-short.pl";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sample       = $sample_files[0];

    for ( my $start = 0 ; $start < $maximum_size ; $start += $bin_size ) {
      my $end       = $start + $bin_size - 1;
      my $curname   = $sample_name . "_" . $start . "_" . $end;
      my $curresult = $curname . ".blastn.tsv";

      my $pbs_file = $self->get_pbs_filename( $pbs_dir, $curname );
      my $pbs_name = basename($pbs_file);
      my $log      = $self->get_log_filename( $log_dir, $curname );

      my $log_desc = $cluster->get_log_description($log);

      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $curresult );
      print $pbs "perl $blastn $curresult $sample $start $end";
      $self->close_pbs( $pbs, $pbs_file );
      print $sh "\$MYCMD ./$pbs_name \n";
    }
  }
  print $sh "exit 0\n";
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $maximum_size = get_option( $config, $section, "maximum_size", 100 );
  my $bin_size     = get_option( $config, $section, "bin_size",     10 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my @result_files = ();

    for ( my $start = 0 ; $start < $maximum_size ; $start += $bin_size ) {
      my $end       = $start + $bin_size - 1;
      my $curname   = $sample_name . "_" . $start . "_" . $end;
      my $curresult = $curname . ".blastn.tsv";
      push( @result_files, "${result_dir}/$curresult" );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $maximum_size = get_option( $config, $section, "maximum_size", 100 );
  my $bin_size     = get_option( $config, $section, "bin_size",     10 );
  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my @result_files = ();

    for ( my $start = 0 ; $start < $maximum_size ; $start += $bin_size ) {
      my $end      = $start + $bin_size - 1;
      my $curname  = $sample_name . "_" . $start . "_" . $end;
      my $pbs_file = $self->get_pbs_filename( $pbs_dir, $curname );
      push( @result_files, $pbs_file );
    }
    $result->{$sample_name} = \@result_files;
  }

  return $result;
}

1;
