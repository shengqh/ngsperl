#!/usr/bin/perl
package Pipeline::ChIPSeq;

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

our %EXPORT_TAGS = ( 'all' => [qw(performChIPSeq performChIPSeqTask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  if ( !defined $def->{cluster} ) {
    $def->{cluster} = 'slurm';
  }

  if ( !defined $def->{fastq_remove_N} ) {
    $def->{fastq_remove_N} = 0;
  }

  if ( !defined $def->{sra_to_fastq} ) {
    $def->{sra_to_fastq} = 0;
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
  my $bwa_index      = $def->{bwa_index};

  my $cqstools   = $def->{cqstools}   or die "define cqstools first";
  my $picard_jar = $def->{picard_jar} or die "define picard_jar first";

  my $macs2_option = $def->{macs2_option} or die "define macs2_option first";
  my $qctable      = $def->{qc_table}     or die "define qctable first";
  my $genome       = $def->{genome}       or die "define genome first, such like hg19, check R ChIPQC package";

  my $config = {
    general => {
      task_name => $task,
      cluster   => $cluster
    },
  };

  $config = merge( $config, $def );

  my $source_ref = "files";
  my @individual;
  my @summary;

  if ($sra_to_fastq) {
    $config->{sra2fastq} = {
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
      option     => "-n -z",
      extension  => "_trim.fastq.gz",
      source_ref => "files",
      cqstools   => $cqstools,
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
    fastqc_raw => {
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
    fastqc_raw_summary => {
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
    bwa => {
      class      => "Alignment::BWA",
      perform    => 1,
      target_dir => "${target_dir}/bwa",
      option     => "",
      bwa_index  => $bwa_index,
      source_ref => "$source_ref",
      picard_jar => $picard_jar,
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    macs2callpeak => {
      class        => "Chipseq::MACS2Callpeak",
      perform      => 1,
      target_dir   => "${target_dir}/macs2callpeak",
      option       => $macs2_option,
      source_ref   => "bwa",
      groups_ref   => "groups",
      controls_ref => "inputs",
      sh_direct    => 0,
      pbs          => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    chipqc => {
      class         => "QC::ChipseqQC",
      perform       => 1,
      target_dir    => "${target_dir}/chipqc",
      option        => "",
      source_ref    => "bwa",
      qctable       => $qctable,
      peaks_ref     => [ "macs2callpeak", "_summits.bed\$" ],
      peak_software => "bed",
      genome        => $genome,
      sh_direct     => 0,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    }
  };

  $config = merge( $config, $processing );
  push @individual, ( "fastqc_raw", "bwa" );
  push @summary, ( "fastqc_raw_summary", "macs2callpeak", "chipqc" );

  $config->{sequencetask} = {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => \@individual,
      step2 => \@summary,
    },
    sh_direct => 0,
    cluster   => $cluster,
    pbs       => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=" . $def->{max_thread},
      "walltime" => $def->{sequencetask_run_time},
      "mem"      => "40gb"
    },
  };

  return ($config);
}

sub performChIPSeq {
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

sub performChIPSeqTask {
  my ( $def, $task ) = @_;

  my $config = getConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
