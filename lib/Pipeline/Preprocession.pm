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

  initDefaultValue( $def, "cluster",                   'slurm' );
  initDefaultValue( $def, "sra_to_fastq",              0 );
  initDefaultValue( $def, "merge_fastq",               0 );
  initDefaultValue( $def, "fastq_remove_N",            0 );
  initDefaultValue( $def, "perform_fastqc",            1 );
  initDefaultValue( $def, "perform_cutadapt",          0 );
  initDefaultValue( $def, "remove_sequences",          "" );
  initDefaultValue( $def, "table_vis_group_text_size", '10' );
  initDefaultValue( $def, "max_thread",                '8' );
  initDefaultValue( $def, "sequencetask_run_time",     '12' );

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
  my $is_paired = $def->{is_paired};
  if ( not defined $is_paired ) {
    $is_paired = $def->{pairend};
  }

  #general
  my $cluster  = getValue( $def, "cluster" );
  my $task     = getValue( $def, "task_name" );
  my $email    = getValue( $def, "email" );
  my $cqstools = getValue( $def, "cqstools" );

  #task
  if ( $def->{sra_to_fastq} ) {
    defined $is_paired or die "Define is_paired first!";
  }

  if ( $def->{merge_fastq} ) {
    defined $is_paired or die "Define is_paired first!";
  }

  my $fastq_remove_N   = getValue( $def, "fastq_remove_N" );
  my $remove_sequences = getValue( $def, "remove_sequences" );    #remove contamination sequences from sequence kit before adapter trimming
  my $run_cutadapt     = getValue( $def, "perform_cutadapt" );
  if ($run_cutadapt) {
    if ( not defined $def->{cutadapt} ) {
      defined $def->{adapter_5} or defined $def->{adapter_3} or getValue( $def, "adapter" );
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
  }

  #data
  my $config = {
    general => {
      task_name => $task,
      cluster   => $cluster
    },
    groups => $def->{groups},
    pairs  => $def->{pairs},
  };

  my $source_ref;
  my $individual = [];
  my $summary    = [];

  if ( $def->{sra_to_fastq} ) {
    $config->{sra2fastq} = {
      class      => "SRA::FastqDump",
      perform    => 1,
      ispaired   => $is_paired,
      target_dir => $def->{target_dir} . "/" . getNextFolderIndex($def) . "sra2fastq",
      option     => "",
      source     => $def->{files},
      sh_direct  => 1,
      cluster    => $def->{cluster},
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    };
    $source_ref = "sra2fastq";
    push @$individual, ("sra2fastq");
  }
  else {
    $config->{files} = getValue( $def, "files" );
    $source_ref = "files";
  }

  if ( $def->{merge_fastq} ) {
    $config->{merge_fastq} = {
      class      => "Format::MergeFastq",
      perform    => 1,
      target_dir => $def->{target_dir} . "/" . getNextFolderIndex($def) . "merge_fastq",
      option     => "",
      source_ref => $source_ref,
      sh_direct  => 1,
      is_paired  => $is_paired,
      cluster    => $def->{cluster},
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=2",
        "walltime" => "2",
        "mem"      => "10gb"
      }
    };
    $source_ref = "merge_fastq";
    push @$individual, ("merge_fastq");
  }

  if ($fastq_remove_N) {
    $config->{fastq_remove_N} = {
      class      => "CQS::FastqTrimmer",
      perform    => 1,
      target_dir => $def->{target_dir} . "/" . getNextFolderIndex($def) . "fastq_remove_N",
      option     => "",
      extension  => "_trim.fastq.gz",
      source_ref => $source_ref,
      sh_direct  => 1,
      cluster    => $def->{cluster},
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      }
    };
    $source_ref = "fastq_remove_N";
    push @$individual, ("fastq_remove_N");
  }

  if ( $def->{perform_fastqc} ) {
    addFastQC( $config, $def, $individual, $summary, "fastqc_raw", $source_ref, $preprocessing_dir );
  }

  if ( length($remove_sequences) ) {
    $config->{"remove_contamination_sequences"} = {
      class      => "CQS::Perl",
      perform    => 1,
      target_dir => $preprocessing_dir . "/" . getNextFolderIndex($def) . "remove_contamination_sequences",
      option     => $remove_sequences,
      output_ext => "_removeSeq.fastq.gz",
      perlFile   => "removeSequenceInFastq.pl",
      source_ref => $source_ref,
      sh_direct  => 1,
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

    if ( $def->{perform_fastqc} ) {
      addFastQC( $config, $def, $individual, $summary, "fastqc_post_remove", $source_ref, $preprocessing_dir );
    }
  }
  
  my $untrimed_ref = $source_ref;

  if ($run_cutadapt) {
    my $cutadapt_section = {
      "cutadapt" => {
        class                            => "Trimmer::Cutadapt",
        perform                          => 1,
        target_dir                       => $preprocessing_dir . "/" . getNextFolderIndex($def) . "cutadapt",
        option                           => $def->{cutadapt_option},
        source_ref                       => $source_ref,
        adapter                          => $def->{adapter},
        adapter_5                        => $def->{adapter_5},
        adapter_3                        => $def->{adapter_3},
        random_bases_remove_after_trim   => $def->{"fastq_remove_random"},
        random_bases_remove_after_trim_5 => $def->{"fastq_remove_random_5"},
        random_bases_remove_after_trim_3 => $def->{"fastq_remove_random_3"},
        extension                        => "_clipped.fastq",
        is_paired                        => $is_paired,
        sh_direct                        => 0,
        cluster                          => $cluster,
        pbs                              => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "20gb"
        },
      }
    };
    if ( defined $def->{cutadapt} ) {
      $cutadapt_section->{cutadapt} = merge( $def->{cutadapt}, $cutadapt_section->{cutadapt} );
    }
    
    my $fastq_len_dir = $preprocessing_dir . "/" . getNextFolderIndex($def) . "fastq_len";
    my $cutadapt = merge(
      $cutadapt_section,
      {
        "fastq_len" => {
          class      => "CQS::FastqLen",
          perform    => 1,
          target_dir => $fastq_len_dir,
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
          target_dir               => $fastq_len_dir,
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
        }
      }
    );

    $config = merge( $config, $cutadapt );

    push @$individual, ( "cutadapt", "fastq_len" );
    push @$summary, ("fastq_len_vis");

    if ( $def->{perform_fastqc} ) {
      addFastQC( $config, $def, $individual, $summary, "fastqc_post_trim", [ "cutadapt", ".fastq.gz" ], $preprocessing_dir );
    }
    $source_ref = [ "cutadapt", ".fastq.gz" ];
  }

  if ( $def->{perform_fastqc} ) {
    my $fastqc_count_vis_files = undef;
    if ( length($remove_sequences) && $run_cutadapt ) {
      $fastqc_count_vis_files = {
        target_dir         => $config->{fastqc_post_trim}->{target_dir},
        parameterFile2_ref => [ "fastqc_post_remove_summary", ".FastQC.summary.reads.tsv\$" ],
        parameterFile3_ref => [ "fastqc_post_trim_summary", ".FastQC.summary.reads.tsv\$" ],
      };
    }
    elsif ( length($remove_sequences) ) {
      $fastqc_count_vis_files = {
        target_dir         => $config->{fastqc_post_remove}->{target_dir},
        parameterFile2_ref => [ "fastqc_post_remove_summary", ".FastQC.summary.reads.tsv\$" ],
      };
    }
    elsif ($run_cutadapt) {
      $fastqc_count_vis_files = {
        target_dir         => $config->{fastqc_post_trim}->{target_dir},
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
  }

  return ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref );
}

1;
