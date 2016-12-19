#!/usr/bin/perl
package Pipeline::ChipSeqBowtie;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performChipSeqBowtie performChipSeqBowtieTask)] );

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

  my $sra_to_fastq   = $def->{sra_to_fastq};
  my $fastq_remove_N = $def->{fastq_remove_N};
  my $email          = $def->{email};
  my $cqstools       = $def->{cqstools} or die "Define cqstools at definition first";
  my $fasta_file     = $def->{fasta_file} or die "Define fasta_file at definition first";
  my $bowtie_index   = $def->{bowtie_index} or die "Define bowtie_index at definition first";
  my $treatments     = $def->{treatments} or die "Define treatments at definition first";
  my $controls       = $def->{controls} or die "Define controls at definition first";

  my $config = {
    general => {
      task_name => $task,
      cluster   => $cluster
    },
    files      => $def->{files},
    treatments => $treatments,
    controls   => $controls
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
    "bowtie1" => {
      class         => "Alignment::Bowtie1",
      perform       => 1,
      target_dir    => "${target_dir}/bowtie1",
      option        => "-v 1 -m 1 --best --strata",
      fasta_file    => $fasta_file,
      source_ref    => [ "cutadapt", ".fastq.gz\$" ],
      bowtie1_index => $bowtie_index,
      sh_direct     => 0,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    "macs1callpeak" => {
      class        => "Chipseq::MACS",
      perform      => 1,
      target_dir   => "${target_dir}/macs1callpeak",
      option       => "-p 1e-9 -w -S --space=50",
      source_ref   => [ "bowtie1", ".bam\$" ],
      groups_ref   => "treatments",
      controls_ref => "controls",
      sh_direct    => 0,
      pbs          => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    "macs1callpeak_bradner_rose2" => {
      class                => "Chipseq::BradnerRose2",
      perform              => 1,
      target_dir           => "${target_dir}/macs1callpeak_bradner_rose2",
      option               => "",
      source_ref           => "bowtie1",
      groups_ref           => "treatments",
      controls_ref         => "controls",
      pipeline_dir         => "/scratch/cqs/shengq1/local/bin/bradnerlab",
      binding_site_bed_ref => [ "macs1callpeak", ".bed\$" ],
      genome               => "hg19",
      sh_direct            => 1,
      pbs                  => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    }
  };
  push @individual, ( "fastqc_raw", "cutadapt", "fastqc_post_trim", "fastq_len", "bowtie1" );
  push @summary, ( "fastqc_raw_summary", "macs1callpeak", "macs1callpeak_bradner_rose2" );

  $config = merge( $config, $processing );

  my $plot_gff = $def->{plot_gff};
  if ($plot_gff) {

    # "-g HG19 -y uniform -r"
    my $bamplot_option = $def->{bamplot_option} or die "Define bamplot_option at definition first";
    my $plotgroups = $def->{plotgroups};
    if ( !defined $plotgroups ) {
      my $files         = $def->{files};
      my @sortedSamples = sort keys %$files;
      $plotgroups = { $task => \@sortedSamples };
    }
    $config->{plotgroups} = $plotgroups;
    $config->{"bamplot"} = {
      class              => "Visualization::Bamplot",
      perform            => 1,
      target_dir         => "${target_dir}/bamplot",
      option             => $bamplot_option,
      source_ref         => "bowtie1",
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

sub performChipSeqBowtie {
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

sub performChipSeqBowtieTask {
  my ( $def, $task ) = @_;
  my $config = getConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
