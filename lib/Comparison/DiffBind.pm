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
use Pipeline::PipelineUtils;

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
  my $homer_annotation_genome = get_option( $config, $section, "homer_annotation_genome", "" );

  my $script = dirname(__FILE__) . "/DiffBind.r";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);
  my $pbs      = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

  my $defaultTissue = $designtable->{Tissue};
  my $defaultFactor = $designtable->{Factor};

  for my $name ( sort keys %$designtable ) {
    if ( $name eq "Tissue" || $name eq "Factor" ) {
      next;
    }

    my $sampleList        = $designtable->{$name};
    my $defaultNameTissue = getValue($sampleList, "Tissue", "Unknown");
    my $defaultNameFactor = getValue($sampleList, "Factor", "Unknown");
    my $comparisons       = getValue($sampleList, "Comparison");

    my $curdir       = create_directory_or_die( $result_dir . "/" . $name );
    my $compFileName = "${name}.comparison.txt";
    my $compfile     = $curdir . "/" . $compFileName;
    open( my $comp, ">$compfile" ) or die "Cannot create $compfile";
    print $comp "Comparison\tGroup1\tGroup2\n";
    for my $comparison (@$comparisons) {
      print $comp $comparison->[0] . "\t" . $comparison->[2] . "\t" . $comparison->[1] . "\n";
    }
    close($comp);

    my $mapFileName = "${name}.config.txt";
    my $mapfile     = $curdir . "/" . $mapFileName;
    open( my $map, ">$mapfile" ) or die "Cannot create $mapfile";
    print $map "SampleID\tTissue\tFactor\tCondition\tReplicate\tbamReads\tPeaks\tPeakCaller\n";
    for my $sampleName ( sort keys %$sampleList ) {
      if ( $sampleName eq "Tissue" || $sampleName eq "Factor" || $sampleName eq "Comparison" ) {
        next;
      }

      my $entryMap = getValue($sampleList, $sampleName);
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
    my $finalPrefix = $name;
    my $finalFile   = $name . "." . $comparisons->[ scalar(@$comparisons) - 1 ]->[0] . ".sig.tsv";
    print $pbs "
cd $curdir

if [ ! -s $finalFile ]; then
  R --vanilla -f $script --args $mapFileName $compFileName $finalPrefix
fi

";
    if ( $homer_annotation_genome ne "" ) {
      for my $comparison (@$comparisons) {
        my $comparisonName = $comparison->[0];
        print $pbs "if [[ -s ${finalPrefix}.${comparisonName}.sig.tsv && ! -s ${finalPrefix}.${comparisonName}.sig.stat.tsv ]]; then 
annotatePeaks.pl ${finalPrefix}.${comparisonName}.sig.tsv $homer_annotation_genome -annStats ${finalPrefix}.${comparisonName}.sig.stat.tsv -go ${finalPrefix}.${comparisonName}.sig.genes.GO > ${finalPrefix}.${comparisonName}.sig.genes.tsv 
fi

";
      }
    }
  }

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );
  my $designtable = get_raw_files( $config, $section, "designtable" );
  my $homer_annotation_genome = get_option( $config, $section, "homer_annotation_genome", "" );

  my $result = {};

  for my $name ( sort keys %$designtable ) {
    my @result_files = ();

    if ( $name eq "Tissue" || $name eq "Factor" ) {
      next;
    }

    my $sampleList        = $designtable->{$name};
    my $comparisons       = $sampleList->{Comparison};

    my $curdir       = create_directory_or_die( $result_dir . "/" . $name );

    my $finalPrefix = $name;
    for my $comparison (@$comparisons){
      my $comparisonName = $comparison->[0];
      my $finalFile   = $name . "." . $comparisonName . ".sig.tsv";
      push (@result_files, $curdir . "/" . $finalFile);
      if ( $homer_annotation_genome ne "" ) {
        push (@result_files, $curdir . "/" . $name . "." . $comparisonName . ".sig.stat.tsv");
        push (@result_files, $curdir . "/" . $name . "." . $comparisonName . ".sig.genes.tsv");
      }
    }
    $result->{$name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
