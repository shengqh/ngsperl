#!/usr/bin/perl
package Pipeline::Preprocession;

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

our %EXPORT_TAGS = ( 'all' => [qw(getPreprocessionConfig)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  initValue( $def, "cluster",                   'slurm' );
  initValue( $def, "sra_to_fastq",              0 );
  initValue( $def, "fastq_remove_N",            0 );
  initValue( $def, "remove_sequences",          "" );
  initValue( $def, "table_vis_group_text_size", '10' );
  initValue( $def, "max_thread",                '8' );
  initValue( $def, "sequencetask_run_time",     '12' );

  return $def;
}

sub getPreprocessionConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = create_directory_or_die( getValue( $def, "target_dir" ) );

  my $preprocessing_dir = $target_dir;
  if ( $def->{subdir} ) {
    $preprocessing_dir = create_directory_or_die( $target_dir . "/preprocessing" );
  }

  $def = initializeDefaultOptions($def);

  #general
  my $cluster  = getValue( $def, "cluster" );
  my $task     = getValue( $def, "task_name" );
  my $email    = getValue( $def, "email" );
  my $cqstools = getValue( $def, "cqstools" );

  #task
  my $sra_to_fastq     = getValue( $def, "sra_to_fastq" );
  my $fastq_remove_N   = getValue( $def, "fastq_remove_N" );
  my $remove_sequences = getValue( $def, "remove_sequences" );    #remove contamination sequences from sequence kit before adapter trimming
  my $run_cutadapt     = getValue( $def, "perform_cutadapt" );
  if ($run_cutadapt) {
    getValue( $def, "adapter" );
    getValue( $def, "min_read_length" );
    my $cutadapt_option = getValue( $def, "cutadapt_option", "" );
    if ( $fastq_remove_N && ( $cutadapt_option !~ /trim\-n/ ) ) {
      $cutadapt_option = $cutadapt_option . " --trim-n";
    }
    $fastq_remove_N = 0;
    if ( $cutadapt_option !~ /\-m/ ) {
      $cutadapt_option = $cutadapt_option . " -m " . $def->{min_read_length};
    }
    $def->{cutadapt_option} = $cutadapt_option;
  }

  #data
  my $files = getValue( $def, "files" );

  my $config = {
    general => {
      task_name => $task,
      cluster   => $cluster
    },
    files => $files,
  };

  my $source_ref = "files";
  my $individual = [];
  my $summary    = [];

  if ($sra_to_fastq) {
    $config->{"sra2fastq"} = {
      class      => "SRA::FastqDump",
      perform    => 1,
      ispaired   => 0,
      target_dir => $preprocessing_dir . "/sra2fastq",
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
    push @$individual, ("sra2fastq");
  }

  if ($fastq_remove_N) {
    $config->{fastq_remove_N} = {
      class      => "CQS::FastqTrimmer",
      perform    => 1,
      target_dir => $preprocessing_dir . "/fastq_remove_N",
      option     => "-n -z",
      extension  => "_trim.fastq.gz",
      source_ref => $source_ref,
      cqstools   => $cqstools,
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
    push @$individual, "fastq_remove_N";
  }

  addFastQC( $config, $def, $individual, $summary, "fastqc_raw", $source_ref, $preprocessing_dir );

  if ( length($remove_sequences) ) {
    $config->{"remove_contamination_sequences"} = {
      class      => "CQS::Perl",
      perform    => 1,
      target_dir => $preprocessing_dir . "/remove_contamination_sequences",
      option     => $remove_sequences,
      output_ext => "_removeSeq.fastq.gz",
      perlFile   => "removeSequenceInFastq.pl",
      source_ref => $source_ref,
      sh_direct  => 0,
      cluster    => $cluster,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "20gb"
      },
    };
    push @$individual, ("remove_contamination_sequences");
    $source_ref = [ "remove_contamination_sequences", ".fastq.gz" ];

    addFastQC( $config, $def, $individual, $summary, "fastqc_post_remove", $source_ref, $preprocessing_dir );
  }

  if ($run_cutadapt) {

    my $cutadapt = {
      "cutadapt" => {
        class                          => "Trimmer::Cutadapt",
        perform                        => 1,
        target_dir                     => $preprocessing_dir . "/cutadapt",
        option                         => $def->{cutadapt_option},
        source_ref                     => $source_ref,
        adapter                        => $def->{adapter},
        extension                      => "_clipped.fastq",
        random_bases_remove_after_trim => $def->{"fastq_remove_random"},
        pairend                        => $def->{pairend},
        sh_direct                      => 0,
        cluster                        => $cluster,
        pbs                            => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "20gb"
        },
      },
      "fastq_len" => {
        class      => "CQS::FastqLen",
        perform    => 1,
        target_dir => $preprocessing_dir . "/fastq_len",
        option     => "",
        source_ref => "cutadapt",
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
      "fastq_len_vis" => {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $preprocessing_dir . "/fastq_len",
        rtemplate                => "countTableVisFunctions.R,fastqLengthVis.R",
        output_file              => ".lengthDistribution",
        output_file_ext          => ".csv",
        parameterSampleFile1_ref => [ "fastq_len", ".len\$" ],
        parameterSampleFile2     => $def->{groups},
        sh_direct                => 1,
        pbs                      => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },
    };

    $config = merge( $config, $cutadapt );

    push @$individual, ( "cutadapt", "fastq_len" );
    push @$summary, ("fastq_len_vis");

    addFastQC( $config, $def, $individual, $summary, "fastqc_post_trim", [ "cutadapt", ".fastq.gz" ], $preprocessing_dir );
    $source_ref = [ "cutadapt", ".fastq.gz" ];
  }

  my $fastqc_count_vis_files = undef;
  if ( length($remove_sequences) && $run_cutadapt ) {
    $fastqc_count_vis_files = {
      target_dir         => $preprocessing_dir . "/fastqc_post_trim",
      parameterFile2_ref => [ "fastqc_post_remove_summary", ".FastQC.summary.reads.tsv\$" ],
      parameterFile3_ref => [ "fastqc_post_trim_summary", ".FastQC.summary.reads.tsv\$" ],
    };
  }
  elsif ( length($remove_sequences) ) {
    $fastqc_count_vis_files = {
      target_dir         => $preprocessing_dir . "/fastqc_post_remove",
      parameterFile2_ref => [ "fastqc_post_remove_summary", ".FastQC.summary.reads.tsv\$" ],
    };
  }
  elsif ($run_cutadapt) {
    $fastqc_count_vis_files = {
      target_dir         => $preprocessing_dir . "/fastqc_post_trim",
      parameterFile2_ref => [ "fastqc_post_trim_summary", ".FastQC.summary.reads.tsv\$" ],
    };
  }

  if ( defined $fastqc_count_vis_files ) {
    $config->{"fastqc_count_vis"} = merge(
      {
        class              => "CQS::UniqueR",
        perform            => 1,
        rtemplate          => "countInFastQcVis.R",
        output_file        => ".countInFastQcVis.Result",
        output_file_ext    => ".Reads.csv",
        sh_direct          => 1,
        parameterFile1_ref => [ "fastqc_raw_summary", ".FastQC.summary.reads.tsv\$" ],
        pbs                => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },
      $fastqc_count_vis_files
    );
    push @$summary, ("fastqc_count_vis");
  }

  #print Dumper($config);

  return ( $config, $individual, $summary, $source_ref, $preprocessing_dir );
}

1;
