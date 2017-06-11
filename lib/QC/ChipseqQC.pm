#!/usr/bin/perl
package QC::ChipseqQC;

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
use Pipeline::PipelineUtils;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_cqc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $bamfiles   = get_raw_files( $config, $section );
  my $qctable    = get_raw_files( $config, $section, "qctable" );
  my $peaksfiles = get_raw_files( $config, $section, "peaks" );
  my $peakSoftware = get_option( $config, $section, "peak_software" );
  my $genome       = get_option( $config, $section, "genome" );
  my $combined     = get_option( $config, $section, "combined" );

  my $script = dirname(__FILE__) . "/ChipseqQC.r";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);
  
  my $expectFiles = $self->result($config, $section);
  my @sortedKeys = (sort keys %$expectFiles);
  my $final_file = $expectFiles->{$sortedKeys[scalar(@sortedKeys)-1]}->[0];
  
  my $pbs      = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  my $mapFiles = writeDesignTable( $result_dir, $section, $qctable, $bamfiles, $peaksfiles, $peakSoftware, $combined, $task_name );
  
  if ($combined) {
    my $mapFileName = $mapFiles->{$task_name};
    print $pbs "R --vanilla -f $script --args $mapFileName $genome \n";
  }
  else {
    for my $qcname ( sort keys %$mapFiles ) {
      my $mapFileName = $mapFiles->{$qcname};
      my $curdir      = $result_dir . "/" . $qcname;
      print $pbs "cd $curdir
R --vanilla -f $script --args $mapFileName $genome \n";
    }
  }

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );
  my $combined = get_option( $config, $section, "combined" );
  my $result = {};

  if ($combined) {
    my @result_files = ();
    my $targetDir    = $result_dir . "/ChIPQCreport";
    push( @result_files, $targetDir . "/ChIPQC.html" );
    $result->{$task_name} = filter_array( \@result_files, $pattern );
  }
  else {
    my $qctable = get_raw_files( $config, $section, "qctable" );
    for my $qcname ( sort keys %$qctable ) {
      if ( $qcname eq "Tissue" || $qcname eq "Factor" ) {
        next;
      }
      my @result_files = ();
      my $curdir       = $result_dir . "/" . $qcname;
      my $targetDir    = $curdir . "/ChIPQCreport";
      push( @result_files, $targetDir . "/ChIPQC.html" );
      $result->{$qcname} = filter_array( \@result_files, $pattern );
    }
  }
  return $result;
}

1;
