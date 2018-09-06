#!/usr/bin/perl
package Pipeline::Parclip;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performParclip performParclipTask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  initDefaultValue( $def, "subdir",                0 );
  initDefaultValue( $def, "sra_to_fastq",          0 );
  initDefaultValue( $def, "bowtie1_option",        "-v 2 -m 10 --best --strata" );
  initDefaultValue( $def, "output_to_same_folder", 1 );
  initDefaultValue( $def, "mappedonly",            1 );

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $email    = getValue( $def, "email" );
  my $cqstools = getValue( $def, "cqstools" );

  $config->{"bowtie1"} = {
    class                 => "Alignment::Bowtie1",
    perform               => 1,
    target_dir            => "${target_dir}/bowtie1",
    option                => getValue( $def, "bowtie1_option" ),
    bowtie1_index         => getValue( $def, "bowtie1_index" ),
    output_to_same_folder => getValue( $def, "output_to_same_folder" ),
    mappedonly            => getValue( $def, "mappedonly" ),
    source_ref            => $source_ref,
    sh_direct             => 1,
    pbs                   => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };
  $config->{"bowtie1_summary"} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/bowtie1",
    rtemplate                => "../Alignment/Bowtie1Summary.r",
    output_file              => "",
    output_file_ext          => ".csv",
    parameterSampleFile1_ref => [ "bowtie1", ".log\$" ],
    sh_direct                => 1,
    pbs                      => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "10gb"
    },
  };

  push @$individual, "bowtie1";
  push @$summary, "bowtie1_summary";

  $config->{"PARalyzer"} = {
    class          => "ParClip::PARalyzer",
    perform        => 1,
    target_dir     => "${target_dir}/paralyzer",
    option         => "",
    source_ref     => "bowtie1",
    genome2bit     => getValue( $def, "genome_2bit" ),
    mirna_db       => getValue( $def, "mirna_db" ),
    sorted_by_name => 0,
    sh_direct      => 1,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  };
  push @$individual, "PARalyzer";

  $config->{"sequencetask"} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_1 => $individual,
      step_2 => $summary,
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  return ($config);
}

sub performParclip {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

sub performParclipTask {
  my ( $def, $task ) = @_;
  my $config = getConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
