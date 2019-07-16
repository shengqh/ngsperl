#!/usr/bin/perl
package Pipeline::TGIRTSmallRNA;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::SmallRNAUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(getTGIRTSmallRNAConfig performTGIRTSmallRNA performTGIRTSmallRNATask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub getTGIRTSmallRNAConfig {
  my ($def) = @_;

  my ( $config, $individual_ref, $summary_ref, $cluster, $source_ref, $preprocessing_dir, $class_independent_dir ) = getPrepareConfig( $def, 0 );
  my @individual = @{$individual_ref};
  my @summary    = @{$summary_ref};

  my $file_ref = "files";
  if ( defined $config->{fastq_remove_N} ) {
    $file_ref = "fastq_remove_N";
  }

  my $tgirt = {
    check_cca => {
      class              => "SmallRNA::TGIRTCheckCCA",
      perform            => 1,
      target_dir         => $def->{target_dir} . "/check_cca",
      option             => "",
      source_ref         => [ 'identical', '.fastq.gz$' ],
      untrimmedFastq_ref => $file_ref,
      sh_direct          => 0,
      pbs                => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb"
      },
    },
    fastq_trna => {
      class        => "SmallRNA::TGIRTNTA",
      perform      => 1,
      target_dir   => $def->{target_dir} . "/fastq_trna",
      option       => "",
      extension    => '_clipped_identical_trna.fastq.gz',
      source_ref   => [ 'identical', '.fastq.gz$' ],
      ccaFile_ref  => "check_cca",
      seqcount_ref => [ 'identical', '.dupcount$' ],
      sh_direct    => 0,
      pbs          => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb"
      },
    },
    star_tRNA => {
      class                     => "Alignment::STAR",
      perform                   => 1,
      target_dir                => $def->{target_dir} . "/star_tRNA",
      option                    => "--outSAMunmapped Within --alignIntronMax 117 --outFilterMismatchNoverLmax 0.05 --clip5pNbases 3 --alignEndsType EndToEnd --outSAMattributes NH HI NM MD AS XS",
      source_ref                => [ 'fastq_trna', '.*.gz$' ],
      genome_dir                => $def->{star_index_directory},
      output_sort_by_coordinate => 1,
      output_unsorted           => 1,
      sh_direct                 => 0,
      pbs                       => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "30gb"
      },
    },
    star_otherSmallRNA => {
      class                     => "Alignment::STAR",
      perform                   => 1,
      target_dir                => $def->{target_dir} . "/star_otherSmallRNA",
      option                    => "--outSAMunmapped Within --alignIntronMax 1 --outFilterMismatchNoverLmax 0.05 --alignEndsType EndToEnd --outSAMattributes NH HI NM MD AS XS",
      source_ref                => [ 'identical', '.fastq.gz$' ],
      genome_dir                => $def->{star_index_directory},
      output_sort_by_coordinate => 1,
      output_unsorted           => 1,
      sh_direct                 => 1,
      pbs                       => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "30gb"
      },
    },
    star_tgirt_count => {
      class              => "SmallRNA::TGIRTCount",
      perform            => 1,
      target_dir         => $def->{target_dir} . "/star_tgirt_count",
      option             => "",
      source_ref         => [ "star_tRNA", "_Aligned.out.bam" ],
      other_smallrna_ref => [ "star_otherSmallRNA", "_Aligned.out.bam" ],
      fastq_files_ref    => "identical",
      seqcount_ref       => [ "identical", ".dupcount\$" ],
      coordinate_file    => $def->{coordinate},
      fasta_file         => $def->{coordinate_fasta},
      sh_direct          => 1,
      pbs                => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "30gb"
      },
    },
    star_tgirt_table => {
      class      => "CQS::SmallRNATable",
      perform    => 1,
      target_dir => $def->{target_dir} . "/star_tgirt_table",
      option     => "",
      source_ref => [ "star_tgirt_count", ".mapped.xml" ],
      prefix     => "tigrt_",
      sh_direct  => 1,
      cluster    => $cluster,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },

    star_tgirt_category => {
      class      => 'CQS::SmallRNACategory',
      perform    => 1,
      target_dir => $def->{target_dir} . "/star_tgirt_category",
      option     => '',
      source_ref => [ 'star_tgirt_count', '.info$' ],
      sh_direct  => 1,
      pbs        => {
        email    => $def->{email},
        walltime => '72',
        mem      => '40gb',
        nodes    => '1:ppn=1'
      },
    },
  };

  push @individual, ( "check_cca", "fastq_trna", "star_tRNA", "star_otherSmallRNA", "star_tgirt_count" );
  push @summary, ( "star_tgirt_table", "star_tgirt_category" );

  $config = merge( $config, $tgirt );
  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => $def->{target_dir} . '/sequencetask',
    option     => '',
    source     => {
      step1 => \@individual,
      step2 => \@summary,
    },
    sh_direct => 0,
    cluster   => 'slurm',
    pbs       => {
      'email'    => $def->{email},
      'nodes'    => '1:ppn=' . $def->{max_thread},
      'walltime' => '72',
      'mem'      => '40gb'
    },
  };

  return ($config);
}

sub performTGIRTSmallRNA {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getTGIRTSmallRNAConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

sub performTGIRTSmallRNATask {
  my ( $def, $task ) = @_;

  my $config = getTGIRTSmallRNAConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
