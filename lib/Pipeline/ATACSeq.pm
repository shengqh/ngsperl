#!/usr/bin/perl
package Pipeline::ATACSeq;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::Preprocession;
use Pipeline::PipelineUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performATACSeq performATACSeqTask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  initValue( $def, "cluster",               "slurm" );
  initValue( $def, "max_thread",            8 );
  initValue( $def, "sequencetask_run_time", 12 );

  #Tasks
  initValue( $def, "sra_to_fastq", 0 );

  initValue( $def, "fastq_remove_N", 0 );

  initValue( $def, "perform_cutadapt", 0 );
  if ( getValue( $def, "perform_cutadapt" ) ) {
    initValue( $def, "adapter",         "CTGTCTCTTATACACATCT" );
    initValue( $def, "min_read_length", 36 );
    initValue( $def, "cutadapt_option", "-q 30" );
  }

  initValue( $def, "perform_rose",          0 );
  initValue( $def, "perform_coltron",       0 );
  initValue( $def, "annotate_nearest_gene", 0 );

  initValue( $def, "perform_homer_motifs", 0 );
  if ( getValue( $def, "perform_homer_motifs" ) ) {
    initValue( $def, "homer_option", "" );
    initValue( $def, "homer_genome", "hg19" );
  }

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir ) = getPreprocessionConfig($def);

  my $cluster = $def->{cluster};
  my $task    = $def->{task_name};

  my $email = getValue( $def, "email" );
  my $macs2call_option = getValue( $def, "macs2call_option", "-f BEDPE --broad -g hs -B -q 0.01 --broad-cutoff 0.01 --nomodel --slocal 20000 --llocal 20000 --keep-dup all" );

  my $perform_rose = getValue( $def, "perform_rose" );
  my $perform_coltron = 0;
  if ($perform_rose) {
    $perform_coltron = getValue( $def, "perform_coltron" );
  }

  my $gene_bed;
  my $annotate_nearest_gene = getValue( $def, "annotate_nearest_gene" );
  if ($annotate_nearest_gene) {
    $gene_bed = getValue( $def, "gene_bed" );
  }

  $config->{treatments} = getValue( $def, "treatments" );

  my $processing = {
    "bwa" => {
      class              => "Alignment::BWA",
      perform            => 1,
      target_dir         => "${target_dir}/bwa",
      option             => "",
      bwa_index          => getValue( $def, "bwa_fasta" ),
      picard_jar         => getValue( $def, "picard_jar" ),
      source_ref         => $source_ref,
      sort_by_coordinate => 1,
      sh_direct          => 0,
      pbs                => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    "bwa_cleanbam" => {
      class                   => "ATACseq::CleanBam",
      perform                 => 1,
      target_dir              => "${target_dir}/bwa_cleanbam",
      option                  => "",
      source_ref              => "bwa",
      picard_jar              => getValue( $def, "picard_jar" ),
      remove_chromosome       => "M",
      keep_chromosome         => "chr",
      is_sorted_by_coordinate => 1,
      sh_direct               => 0,
      pbs                     => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "240",
        "mem"      => "40gb"
      },
    },
    "bwa_bam2bed" => {
      class                   => "Format::Bam2Bed",
      perform                 => 1,
      target_dir              => "${target_dir}/bwa_bam2bed",
      option                  => "",
      source_ref              => "bwa_cleanbam",
      blacklist_file          => "/scratch/cqs/shengq1/references/mappable_region/hg19/wgEncodeDacMapabilityConsensusExcludable.bed",
      is_sorted_by_name       => 0,
      is_paired_end           => 1,
      maximum_fragment_length => 1000,
      minimum_fragment_length => 30,
      sh_direct               => 1,
      pbs                     => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
  };
  push @$individual, ( "bwa", "bwa_cleanbam", "bwa_bam2bed" );
  if ( $perform_rose || ( !defined $def->{replicates} ) ) {
    $config->{"bwa_macs2callpeak"} = {
      class      => "Chipseq::MACS2Callpeak",
      perform    => 1,
      target_dir => "${target_dir}/bwa_macs2callpeak",
      option     => $macs2call_option,
      source_ref => "bwa_bam2bed",
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    push @$summary, ("bwa_macs2callpeak");

    if ( getValue( $def, "perform_homer_motifs" ) ) {
      $config->{bwa_macs2callpeak_homer_motifs} = {
        class        => "Homer::FindMotifs",
        option       => getValue( $def, "homer_option" ),
        perform      => 1,
        source_ref   => [ "bwa_macs2callpeak", ".bed\$" ],
        target_dir   => $target_dir . "/bwa_macs2callpeak_homer_motifs",
        homer_genome => getValue( $def, "homer_genome" ),
        sh_direct    => 1,
        pbs          => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push @$summary, ("bwa_macs2callpeak_homer_motifs");
    }
  }

  if ($perform_rose) {
    $config->{"bwa_macs2callpeak_bradner_rose"} = {
      class                => "Chipseq::BradnerRose2",
      perform              => 1,
      target_dir           => "${target_dir}/bwa_macs2callpeak_bradner_rose",
      option               => "",
      source_ref           => "bwa_cleanbam",
      groups_ref           => "treatments",
      pipeline_dir         => "/scratch/cqs/shengq1/local/bin/bradnerlab",
      binding_site_bed_ref => [ "bwa_macs2callpeak", ".bed\$" ],
      genome               => "hg19",
      sh_direct            => 1,
      pbs                  => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    push @$summary, ("bwa_macs2callpeak_bradner_rose");

    if ($perform_coltron) {
      $config->{"bwa_macs2callpeak_bradner_rose_coltron"} = {
        class              => "Chipseq::Coltron",
        perform            => 1,
        target_dir         => "${target_dir}/bwa_macs2callpeak_bradner_rose_coltron",
        option             => "",
        source_ref         => "bwa_cleanbam",
        groups_ref         => "treatments",
        enhancer_files_ref => [ "bwa_macs2callpeak_bradner_rose", "_AllEnhancers.table.txt" ],
        genome             => "HG19",
        pipeline_dir       => "/scratch/cqs/shengq1/local/bin/bradnerlab",
        sh_direct          => 1,
        pbs                => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      };
      push @$summary, ("bwa_macs2callpeak_bradner_rose_coltron");
    }
  }

  if ( defined $def->{comparison} ) {
    my $peakTask = "bwa_macs2callpeak";
    if ( defined $def->{replicates} ) {
      $config->{"bwa_macs2callpeak_replicates"} = {
        class      => "Chipseq::MACS2Callpeak",
        perform    => 1,
        target_dir => "${target_dir}/bwa_macs2callpeak_replicates",
        option     => $macs2call_option,
        source_ref => "bwa_bam2bed",
        groups     => $def->{"replicates"},
        sh_direct  => 0,
        pbs        => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      };
      push @$summary, ("bwa_macs2callpeak_replicates");
      $peakTask = "bwa_macs2callpeak_replicates";

      if ( getValue( $def, "perform_homer_motifs" ) ) {
        $config->{bwa_macs2callpeak_replicates_homer_motifs} = {
          class        => "Homer::FindMotifs",
          option       => getValue( $def, "homer_option" ),
          perform      => 1,
          source_ref   => [ "bwa_macs2callpeak_replicates", ".bed\$" ],
          target_dir   => $target_dir . "/bwa_macs2callpeak_replicates_homer_motifs",
          homer_genome => getValue( $def, "homer_genome" ),
          sh_direct    => 1,
          pbs          => {
            "email"     => $def->{email},
            "emailType" => $def->{emailType},
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        };
        push @$summary, ("bwa_macs2callpeak_replicates_homer_motifs");
      }
    }

    $config->{"bwa_macs2diff"} = {
      class      => "Chipseq::MACS2Bdgdiff",
      perform    => 1,
      target_dir => "${target_dir}/bwa_macs2diff",
      option     => "",
      source_ref => $peakTask,
      groups     => $def->{comparison},
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    push @$summary, ("bwa_macs2diff");

    if ($annotate_nearest_gene) {
      $config->{bwa_macs2diff_gene} = {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $target_dir . "/bwa_macs2diff",
        rtemplate                => "../Annotation/findNearestGene.r",
        output_file              => "",
        output_file_ext          => ".Category.Table.csv",
        parameterSampleFile1_ref => "bwa_macs2diff",
        parameterFile1           => $gene_bed,
        rCode                    => '',
        sh_direct                => 1,
        pbs                      => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push @$summary, ("bwa_macs2diff_gene");
    }

    if ( getValue( $def, "perform_homer_motifs" ) ) {
      $config->{bwa_macs2diff_homer_motifs} = {
        class        => "Homer::FindMotifs",
        option       => getValue( $def, "homer_option" ),
        perform      => 1,
        source_ref   => "bwa_macs2diff",
        target_dir   => $target_dir . "/bwa_macs2diff_homer_motifs",
        homer_genome => getValue( $def, "homer_genome" ),
        sh_direct    => 1,
        pbs          => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push @$summary, ("bwa_macs2diff_homer_motifs");
    }
  }

  $config = merge( $config, $processing );

  my $plot_gff = $def->{plot_gff};
  if ($plot_gff) {

    # "-g HG19 -y uniform -r"
    my $bamplot_option = getValue( $def, "bamplot_option" );
    my $plotgroups = $def->{plotgroups};
    if ( !defined $plotgroups ) {
      my $files         = $def->{files};
      my @sortedSamples = sort keys %$files;
      $plotgroups = { $task => \@sortedSamples };
    }
    $config->{"plotgroups"} = $plotgroups;
    $config->{"bamplot"}    = {
      class              => "Visualization::Bamplot",
      perform            => 1,
      target_dir         => "${target_dir}/bamplot",
      option             => $bamplot_option,
      source_ref         => "bwa",
      groups_ref         => "plotgroups",
      gff_file           => $plot_gff,
      is_rainbow_color   => 0,
      is_draw_individual => 0,
      is_single_pdf      => 1,
      sh_direct          => 1,
      pbs                => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "1",
        "mem"      => "10gb"
      },
    };
    push @$summary, ("bamplot");
  }

  $config->{"sequencetask"} = {
    class      => "CQS::SequenceTask",
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

sub performATACSeq {
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

sub performATACSeqTask {
  my ( $def, $task ) = @_;
  my $config = getConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
