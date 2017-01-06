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

  my $defaultTissue = delete $qctable->{Tissue};

  my $script = dirname(__FILE__) . "/ChipseqQC.r";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);
  my $pbs      = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

  for my $qcname ( sort keys %$qctable ) {
    my $sampleList          = $qctable->{$qcname};
    my $defaultQcnameTissue = delete $sampleList->{Tissue};
    my $curdir              = create_directory_or_die( $result_dir . "/" . $qcname );
    my $mapfile             = $curdir . "/${qcname}.txt";
    open( my $map, ">$mapfile" ) or die "Cannot create $mapfile";
    print $map "SampleID\tTissue\tFactor\tReplicate\tbamReads\tControlID\tbamControl\tPeaks\tPeakCaller\n";
    for my $sampleName ( sort keys %$sampleList ) {
      my $entryMap = $sampleList->{$sampleName};
      my $sampleId = $entryMap->{Sample} or die "Define Sample for $sampleName in qctable of section $section";
      my $tissue   = $entryMap->{Tissue};
      if ( !defined $tissue ) {
        $tissue = $defaultQcnameTissue;
        if ( !defined $tissue ) {
          $tissue = $defaultTissue;
          if ( !defined $tissue ) {
            die "Define Tissue in $sampleName, or in $qcname, or in section $section";
          }
        }
      }
      my $factor    = $entryMap->{Factor}    or die "Define Factor for $sampleName in qctable of section $section";
      my $replicate = $entryMap->{Replicate} or die "Define Replicate for $sampleName in qctable of section $section";
      my $bamReads  = $bamfiles->{$sampleId}->[0];
      my $controlId = $entryMap->{Input}     or die "Define Input for $sampleName in qctable of section $section";
      my $bamControl = $bamfiles->{$controlId}->[0];
      my $peakFile   = $peaksfiles->{$sampleName}->[0];

      print $map $sampleId . "\t" . $tissue . "\t" . $factor . "\t" . $replicate . "\t" . $bamReads . "\t" . $controlId . "\t" . $bamControl . "\t" . $peakFile . "\t" . $peakSoftware . "\n";
    }
    close($map);

    print $pbs "R --vanilla -f $script --args $mapfile $genome \n";
  }

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $result       = {};
  my @result_files = ();
  push( @result_files, $result_dir . "/bamReport.html" );
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
