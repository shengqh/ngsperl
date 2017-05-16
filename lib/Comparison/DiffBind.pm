#!/usr/bin/perl
package Comparison::DiffBind;

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
  $self->{_suffix} = "_db";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $bamfiles    = get_raw_files( $config, $section );
  my $designtable = get_raw_files( $config, $section, "designtable" );
  my $peaksfiles  = get_raw_files( $config, $section, "peaks" );
  my $peakSoftware = get_option( $config, $section, "peak_software" );

  my $script = dirname(__FILE__) . "/DiffBind.r";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);
  my $pbs      = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

  my $defaultTissue = delete $designtable->{Tissue};
  my $defaultFactor = delete $designtable->{Factor};

  for my $name ( sort keys %$designtable ) {
    my $sampleList        = $designtable->{$name};
    my $defaultNameTissue = delete $sampleList->{Tissue};
    my $defaultNameFactor = delete $sampleList->{Factor};

    my $curdir      = create_directory_or_die( $result_dir . "/" . $name );
    my $mapFileName = "${name}.txt";
    my $mapfile     = $curdir . "/" . $mapFileName;
    open( my $map, ">$mapfile" ) or die "Cannot create $mapfile";
    print $map "SampleID\tTissue\tFactor\tCondition\tReplicate\tbamReads\tPeaks\tPeakCaller\n";
    for my $sampleName ( sort keys %$sampleList ) {
      my $entryMap = $sampleList->{$sampleName};
      my $tissue   = $entryMap->{Tissue};
      if ( !defined $tissue ) {
        $tissue = $defaultNameTissue;
        if ( !defined $tissue ) {
          $tissue = $defaultTissue;
          if ( !defined $tissue ) {
            die "Define Tissue in $sampleName, or in $name, or in section $section";
          }
        }
      }
      my $factor = defined $entryMap->{Factor} ? $entryMap->{Factor} : defined $defaultFactor ? $defaultFactor : die "Define Factor for $sampleName in designtable of section $section";
      my $condition = $entryMap->{Condition} or die "Define Condition for $sampleName in designtable of section $section";
      my $replicate = $entryMap->{Replicate} or die "Define Replicate for $sampleName in designtable of section $section";
      my $bamReads  = $bamfiles->{$sampleName}->[0];
      my $peakFile  = $peaksfiles->{$sampleName}->[0];

      print $map $sampleName . "\t" . $tissue . "\t" . $factor . "\t" . $condition . "\t" . $replicate . "\t" . $bamReads . "\t" . $peakFile . "\t" . $peakSoftware . "\n";
    }
    close($map);
    my $finalPrefix = $name . ".result";
    my $finalFile   = $name . ".result.tsv";
    print $pbs "
cd $curdir

if [ ! -s $finalFile ]; then
  R --vanilla -f $script --args $mapFileName $finalPrefix
fi
";
  }

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );
  my $designtable = get_raw_files( $config, $section, "designtable" );

  my $result = {};

  for my $name ( sort keys %$designtable ) {
    my @result_files = ();
    my $curdir       = $result_dir . "/" . $name;
    push( @result_files, "${curdir}/${name}.result.tsv" );
    $result->{$name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
