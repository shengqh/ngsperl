#!/usr/bin/perl
package Pipeline::RNASeq;

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

  my $fastq_remove_N = $def->{fastq_remove_N};
  my $email          = $def->{email};
  my $bwa_index      = $def->{bwa_index};

  my $cqstools   = $def->{cqstools}   or die "define cqstools first";
  my $picard_jar = $def->{picard_jar} or die "define picard_jar first";
  my $spp_r      = $def->{spp_r}      or die "define spp_r first";
  
  my $groups = $def->{groups};
  my $inputs = $def->{inputs};

  my $config = {
    general => {
      task_name => $task,
      cluster   => $cluster
    },
    files   => $def->{files},
    samples => $def->{samples},
  };

  if ( defined $def->{input} ) {
    $config->{input} = $def->{input};
  }

  my $source_ref = "files";
  my @individual;
  my @summary;

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
    spp => {
      class      => "Chipseq::SPP",
      perform    => 1,
      target_dir => "${target_dir}/spp",
      option     => "",
      source_ref => "$source_ref",
      groups => $groups,
      inputs => $inputs,
      spp_r      => $spp_r,
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },

    }
  };
  push @individual, ( "fastqc_raw", "bwa", "spp" );
  push @summary, ("fastqc_raw_summary");
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
