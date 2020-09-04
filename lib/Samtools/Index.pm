#!/usr/bin/perl
package Samtools::Index;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::GroupTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_ix";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $isbamsorted = $config->{$section}{isbamsorted};
  if ( !defined($isbamsorted) ) {
    $isbamsorted = 0;
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

    my $bam_file = $sample_files[0];

    my $bamSortedFile;
    if ($isbamsorted) {
      $bamSortedFile = $bam_file;
    }
    else {
      ( $bamSortedFile, my $bamSorted ) = get_sorted_bam($bam_file);
      print $pbs "if [ ! -s $bamSortedFile ]; then
  echo samtools_sort=`date`
  samtools sort $bam_file $bamSorted 
fi
";
    }

    my $bamIndexFile = $bamSortedFile . ".bai";
    print $pbs "if [ ! -s $bamIndexFile ]; then
  echo samtools_index=`date`
  samtools index $bamSortedFile 
fi
";
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $isbamsorted = $config->{$section}{isbamsorted};
  if ( !defined($isbamsorted) ) {
    $isbamsorted = 0;
  }

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };

    my $bam_file = $sample_files[0];

    my $bamSortedFile;
    if ($isbamsorted) {
      $bamSortedFile = $bam_file;
    }
    else {
      ( $bamSortedFile, my $bamSorted ) = get_sorted_bam($bam_file);
    }

    my $bamIndexFile = $bamSortedFile . ".bai";

    my @result_files = ();
    push( @result_files, $bamIndexFile );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }

  return $result;
}

1;
