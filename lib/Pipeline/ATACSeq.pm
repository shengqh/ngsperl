#!/usr/bin/perl
package Pipeline::ATACSeq;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
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

  if ( !defined $def->{cluster} ) {
    $def->{cluster} = 'slurm';
  }

  if ( !defined $def->{sra_to_fastq} ) {
    $def->{sra_to_fastq} = 0;
  }

  if ( !defined $def->{fastq_remove_N} ) {
    $def->{fastq_remove_N} = 0;
  }

  if ( !defined $def->{table_vis_group_text_size} ) {
    $def->{table_vis_group_text_size} = "10";
  }

  if ( !defined $def->{max_thread} ) {
    $def->{max_thread} = "8";
  }
  if ( !defined $def->{sequencetask_run_time} ) {
    $def->{sequencetask_run_time} = "12";
  }

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my $cluster = $def->{cluster};
  my $task    = $def->{task_name};

  my $sra_to_fastq     = getValue( $def, "sra_to_fastq",     0 );
  my $fastq_remove_N   = getValue( $def, "fastq_remove_N",   0 );
  my $email            = getValue( $def, "email" );
  my $cqstools         = getValue( $def, "cqstools" );
  my $picard_jar       = getValue( $def, "picard_jar" );
  my $bwa_fasta        = getValue( $def, "bwa_fasta" );
  my $treatments       = getValue( $def, "treatments" );
  my $pairend          = getValue( $def, "pairend" );
  my $macs2call_option = getValue( $def, "macs2call_option", "-f BEDPE --broad -g hs -B -q 0.01 --broad-cutoff 0.01 --nomodel --slocal 20000 --llocal 20000 --keep-dup all" );

  my $perform_rose = getValue( $def, "perform_rose" );
  my $perform_coltron = 0;
  if ($perform_rose) {
    $perform_coltron = getValue( $def, "perform_coltron" );
  }

  my $config = {
    general => {
      task_name => $task,
      cluster   => $cluster
    },
    files      => $def->{files},
    treatments => $treatments
  };

  my $source_ref = "files";
  my @individual;
  my @summary;

  if ($sra_to_fastq) {
    $config->{"sra2fastq"} = {
      class      => "SRA::FastqDump",
      perform    => 1,
      ispaired   => 0,
      target_dir => "${target_dir}/sra2fastq",
      option     => "",
      source_ref => "files",
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    };
    $source_ref = "sra2fastq";
    push @individual, ("sra2fastq");
  }

  if ($fastq_remove_N) {
    $config->{fastq_remove_N} = {
      class      => "CQS::FastqTrimmer",
      perform    => $fastq_remove_N,
      target_dir => $target_dir . "/fastq_remove_N",
      option     => "",
      extension  => "_trim.fastq.gz",
      source_ref => "files",
      cluster    => $cluster,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      }
    };
    $source_ref = "fastq_remove_N";
    push @individual, ("fastq_remove_N");
  }

  my $processing = {
    "fastqc_raw" => {
      class      => "QC::FastQC",
      perform    => 1,
      target_dir => $target_dir . "/fastqc_raw",
      option     => "",
      source_ref => $source_ref,
      cluster    => $cluster,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      },
    },
    "fastqc_raw_summary" => {
      class      => "QC::FastQCSummary",
      perform    => 1,
      target_dir => $target_dir . "/fastqc_raw",
      cqstools   => $cqstools,
      option     => "",
      cluster    => $cluster,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      },
    },
    "cutadapt" => {
      class      => "Trimmer::Cutadapt",
      perform    => 1,
      target_dir => "${target_dir}/cutadapt",
      option     => "-m 30 --trim-n",
      source_ref => $source_ref,
      adapter    => "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",    #trueseq adapter
      extension  => "_clipped.fastq",
      pairend    => $pairend,
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    "fastqc_post_trim" => {
      class      => "QC::FastQC",
      perform    => 1,
      target_dir => "${target_dir}/fastqc_post_trim",
      option     => "",
      sh_direct  => 1,
      source_ref => [ "cutadapt", ".fastq.gz" ],
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      },
    },
    "fastqc_post_trim_summary" => {
      class      => "QC::FastQCSummary",
      perform    => 1,
      sh_direct  => 1,
      target_dir => "${target_dir}/fastqc_post_trim",
      cqstools   => $cqstools,
      option     => "",
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      },
    },
    "fastq_len" => {
      class      => "CQS::FastqLen",
      perform    => 1,
      target_dir => "$target_dir/fastq_len",
      option     => "",
      source_ref => "cutadapt",
      cqstools   => $cqstools,
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    "bwa" => {
      class              => "Alignment::BWA",
      perform            => 1,
      target_dir         => "${target_dir}/bwa",
      option             => "",
      bwa_index          => $bwa_fasta,
      picard_jar         => $picard_jar,
      source_ref         => [ "cutadapt", ".fastq.gz" ],
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
      picard_jar              => $picard_jar,
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
    "bwa_macs2callpeak" => {
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
    },
  };
  push @individual, ( "fastqc_raw", "cutadapt", "fastqc_post_trim", "fastq_len", "bwa", "bwa_cleanbam", "bwa_bam2bed" );
  push @summary, ( "fastqc_raw_summary", "bwa_macs2callpeak" );

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
    push @summary, ("bwa_macs2callpeak_bradner_rose");

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
      push @summary, ("bwa_macs2callpeak_bradner_rose_coltron");
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
      push @summary, ("bwa_macs2callpeak_replicates");
      $peakTask = "bwa_macs2callpeak_replicates";
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
    push @summary, ("bwa_macs2diff");
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
    push @summary, ("bamplot");
  }

  $config->{"sequencetask"} = {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_1 => \@individual,
      step_2 => \@summary,
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
