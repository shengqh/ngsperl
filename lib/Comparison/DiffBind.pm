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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $bamfiles    = get_raw_files( $config, $section );
  my $designtable = get_raw_files( $config, $section, "design_table" );
  my $peaksfiles  = get_raw_files( $config, $section, "peaks" );
  my $peakSoftware = get_option( $config, $section, "peak_software" );
  my $homer_annotation_genome = get_option( $config, $section, "homer_annotation_genome", "" );

  my $treatments = $config->{$section}{"groups"};
  my $controls = $config->{$section}{"inputs"};
  if ( !defined $controls ) {
    $controls = $config->{$section}{"controls"};
  }

  my $script = dirname(__FILE__) . "/DiffBind.r";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);
  my $pbs      = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

  my $mapFiles = writeDesignTable( $result_dir, $section, $designtable, $bamfiles, $peaksfiles, $peakSoftware, 0, $task_name, $treatments, $controls );

  for my $name ( sort keys %$mapFiles ) {
    my $mapFileName = basename($mapFiles->{$name});

    my $sampleList = $designtable->{$name};
    my $comparisons = getValue( $sampleList, "Comparison" );
    my $minOverlap = getValue( $sampleList, "MinOverlap");

    my $curdir       = $result_dir . "/" . $name;
    my $compFileName = "${name}.comparison.txt";
    my $compfile     = $curdir . "/" . $compFileName;
    open( my $comp, ">$compfile" ) or die "Cannot create $compfile";
    print $comp "Comparison\tGroup1\tGroup2\n";
    for my $comparison (@$comparisons) {
      print $comp $comparison->[0] . "\t" . $comparison->[2] . "\t" . $comparison->[1] . "\n";
    }
    close($comp);

    my $overlapFileName = "${name}.minoverlap.txt";

    my $finalPrefix = $name;
    my $finalFile   = $name . "." . $comparisons->[ scalar(@$comparisons) - 1 ]->[0] . ".sig.tsv";
    print $pbs "
cd $curdir

if [ ! -s $finalFile ]; then
  R --vanilla -f $script --args $mapFileName $compFileName $overlapFileName $finalPrefix
fi

";

    if ( $homer_annotation_genome ne "" ) {
      for my $oKey (sort keys %$minOverlap){
        my $peaksFile = $oKey . ".peaks";
        print $pbs "annotatePeaks.pl <(awk -F'\\t' -v OFS='\\t' '{print \$1,\$2,\$3,\"peak\"\$1\":\"\$2,\"NA\",\"+\"}' $peaksFile) $homer_annotation_genome -annStats ${oKey}.stat.tsv -go ${oKey}.genes.GO > ${oKey}.genes.tsv \n\n";
      }

      for my $comparison (@$comparisons) {
        my $comparisonName = $comparison->[0];
        print $pbs "annotatePeaks.pl ${finalPrefix}.${comparisonName}.sig.tsv $homer_annotation_genome -annStats ${finalPrefix}.${comparisonName}.sig.stat.tsv -go ${finalPrefix}.${comparisonName}.sig.genes.GO > ${finalPrefix}.${comparisonName}.sig.genes.tsv \n\n";
      }
    }
  }

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  my $designtable = get_raw_files( $config, $section, "design_table" );
  my $homer_annotation_genome = get_option( $config, $section, "homer_annotation_genome", "" );

  my $result = {};

  for my $name ( sort keys %$designtable ) {
    if ( $name eq "Tissue" || $name eq "Factor" ) {
      next;
    }

    my @result_files = ();

    my $sampleList  = $designtable->{$name};
    my $comparisons = $sampleList->{Comparison};
    my $minOverlap = $sampleList->{MinOverlap};

    my $curdir = $result_dir . "/" . $name;

    my $finalPrefix = $name;
    for my $comparison (@$comparisons) {
      my $comparisonName = $comparison->[0];
      push( @result_files, "${curdir}/${name}.${comparisonName}.sig.tsv" );
      push( @result_files, "${curdir}/${name}.${comparisonName}.sig.bed" );
      if ( $homer_annotation_genome ne "" ) {
        push( @result_files, "${curdir}/${name}.${comparisonName}.sig.stat.tsv" );
        push( @result_files, "${curdir}/${name}.${comparisonName}.sig.genes.tsv" );
      }
    }
    for my $oKey (sort keys %$minOverlap){
      push( @result_files, $curdir . "/" . $oKey . ".peaks" );
      push( @result_files, $curdir . "/" . $oKey . ".pdf" );
    }
    push( @result_files, $curdir . "/Overlap-condition.pdf" );

    $result->{$name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
