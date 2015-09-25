#!/usr/bin/perl
package Pipeline::SmallRNAUtils;

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

our %EXPORT_TAGS = ( 'all' => [qw(getPrepareConfig)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub getPrepareConfig {
  my ($def) = @_;

  create_directory_or_die( $def->{target_dir} );

  my $cluster = $def->{cluster};
  if ( !defined $cluster ) {
    $cluster = "slurm";
  }

  my $fastq_remove_N = $def->{fastq_remove_N};
  my $run_cutadapt   = $def->{run_cutadapt};

  my $config = {
    general => {
      task_name => $def->{task_name},
      cluster   => $cluster
    },
    files => $def->{files}
  };

  my @individual = ();
  my @summary    = ();

  my $source_ref = "files";
  my $len_ref = "files";
  if ( !defined $fastq_remove_N || $fastq_remove_N ) {
    $config->{fastq_remove_N} = {
      class      => "CQS::FastqTrimmer",
      perform    => $fastq_remove_N,
      target_dir => $def->{target_dir} . "/fastq_remove_N",
      option     => "-n -z",
      extension  => "_trim.fastq.gz",
      source_ref => "files",
      cqstools   => $def->{cqstools},
      cluster    => $cluster,
      sh_direct  => 1,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      }
    };
    $source_ref = "fastq_remove_N";
    $len_ref = "fastq_remove_N";
    push @individual, "fastq_remove_N";
  }

  my $qc = {};
  if ( !defined $run_cutadapt || $run_cutadapt ) {
    my $adapter = $def->{adapter};
    if ( !defined $adapter ) {
      $adapter = "TGGAATTCTCGGGTGCCAAGG";
    }

    $qc = {
      fastqc_pre_trim => {
        class      => "QC::FastQC",
        perform    => 1,
        target_dir => $def->{target_dir} . "/fastqc_pre_trim",
        option     => "",
        source_ref => $source_ref,
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      },
      fastqc_pre_trim_summary => {
        class      => "QC::FastQCSummary",
        perform    => 1,
        target_dir => $def->{target_dir} . "/fastqc_pre_trim",
        cqstools   => $def->{cqstools},
        option     => "",
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      },
      cutadapt => {
        class      => "Cutadapt",
        perform    => 1,
        target_dir => $def->{target_dir} . "/cutadapt",
        option     => "-O 10 -m " . $def->{min_read_length},
        source_ref => $source_ref,
        adapter    => $adapter,
        extension  => "_clipped.fastq",
        sh_direct  => 1,
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "20gb"
        },
      },
      fastqc_post_trim => {
        class      => "QC::FastQC",
        perform    => 1,
        target_dir => $def->{target_dir} . "/fastqc_post_trim",
        option     => "",
        source_ref => [ "cutadapt", ".fastq.gz" ],
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      },
      fastqc_post_trim_summary => {
        class      => "QC::FastQCSummary",
        perform    => 1,
        target_dir => $def->{target_dir} . "/fastqc_post_trim",
        cqstools   => $def->{cqstools},
        option     => "",
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      }
    };
    $source_ref = [ "cutadapt", ".fastq.gz" ];
    $len_ref = "cutadapt";
    push @individual, ( "fastqc_pre_trim", "cutadapt", "fastqc_post_trim" );
    push @summary, ( "fastqc_pre_trim_summary", "fastqc_post_trim_summary" );
  }
  else {
    $qc = {
      fastqc => {
        class      => "QC::FastQC",
        perform    => 1,
        target_dir => $def->{target_dir} . "/fastqc",
        option     => "",
        source_ref => $source_ref,
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      },
      fastqc_summary => {
        class      => "QC::FastQCSummary",
        perform    => 1,
        target_dir => $def->{target_dir} . "/fastqc",
        cqstools   => $def->{cqstools},
        option     => "",
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      }
    };
    push @individual, ("fastqc");
    push @summary,    ("fastqc_summary");
  }

  $config = merge( $config, $qc );

  #print Dumper($config);

  my $preparation = {
    fastq_len => {
      class      => "FastqLen",
      perform    => 1,
      target_dir => $def->{target_dir} . "/fastq_len",
      option     => "",
      source_ref => $len_ref,
      cqstools   => $def->{cqstools},
      sh_direct  => 1,
      cluster    => $cluster,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    identical => {
      class      => "FastqIdentical",
      perform    => 1,
      target_dir => $def->{target_dir} . "/identical",
      option     => "",
      source_ref => $source_ref,
      cqstools   => $def->{cqstools},
      extension  => "_clipped_identical.fastq.gz",
      sh_direct  => 1,
      cluster    => $cluster,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    identical_sequence_count_table => {
      class      => "CQS::SmallRNASequenceCountTable",
      perform    => 1,
      target_dir => $def->{target_dir} . "/identical_sequence_count_table",
      option     => "",
      source_ref => [ "identical", ".dupcount\$" ],
      cqs_tools  => $def->{cqstools},
      suffix     => "_sequence",
      sh_direct  => 1,
      cluster    => $cluster,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    identical_NTA => {
      class        => "CQS::FastqMirna",
      perform      => 1,
      target_dir   => $def->{target_dir} . "/identical_NTA",
      option       => "-l " . $def->{min_read_length},
      source_ref   => [ "identical", ".fastq.gz\$" ],
      seqcount_ref => [ "identical", ".dupcount\$" ],
      cqstools     => $def->{cqstools},
      extension    => "_clipped_identical_NTA.fastq.gz",
      sh_direct    => 1,
      cluster      => $cluster,
      pbs          => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
  };

  push @individual, ( "fastq_len", "identical", "identical_NTA" );
  push @summary, ("identical_sequence_count_table");

  $config = merge( $config, $preparation );

  return ($config, \@individual, \@summary, $cluster);
}

1;
