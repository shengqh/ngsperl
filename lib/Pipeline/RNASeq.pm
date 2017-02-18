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

our %EXPORT_TAGS = ( 'all' => [qw(performRNASeq performRNASeqTask)] );

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

  if ( !defined $def->{perform_rnaseqc} ) {
    $def->{perform_rnaseqc} = 0;
  }

  if ( !defined $def->{perform_qc3bam} ) {
    $def->{perform_qc3bam} = 0;
  }

  if ( !defined $def->{perform_bamplot} ) {
    $def->{perform_bamplot} = 0;
  }

  return $def;
}

sub getRNASeqConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my $cluster = $def->{cluster};
  my $task    = $def->{task_name};

  my $sra_to_fastq               = $def->{sra_to_fastq};
  my $fastq_remove_N             = $def->{fastq_remove_N};
  my $email                      = $def->{email};
  my $cqstools                   = $def->{cqstools} or die "Define cqstools at definition first";
  my $star_index                 = $def->{star_index} or die "Define star_index at definition first";
  my $transcript_gtf             = $def->{transcript_gtf} or die "Define transcript_gtf at definition first";
  my $name_map_file              = $def->{name_map_file} or die "Define tramscript name_map_file at definition first";
  my $use_pearson_in_hca         = defined $def->{use_pearson_in_hca} ? $def->{use_pearson_in_hca} : 1;
  my $top25cv_in_hca             = defined $def->{top25cv_in_hca} ? $def->{top25cv_in_hca} : 0;
  my $use_green_red_color_in_hca = defined $def->{use_green_red_color_in_hca} ? $def->{use_green_red_color_in_hca} : 1;

  if ( $def->{perform_rnaseqc} ) {
    defined $def->{rnaseqc_jar} or die "Define rnaseqc_jar first!";
    ( -e $def->{rnaseqc_jar} )  or die "rnaseqc_jar not exists " . $def->{rnaseqc_jar};
    defined $def->{fasta_file}  or die "Define fasta_file for rnaseqc first!";
    ( -e $def->{fasta_file} )   or die "fasta_file not exists " . $def->{fasta_file};
  }

  if ( $def->{perform_qc3bam} ) {
    defined $def->{qc3_perl} or die "Define qc3_perl first!";
    ( -e $def->{qc3_perl} )  or die "qc3_perl not exists " . $def->{qc3_perl};
  }

  if ( $def->{perform_bamplot} ) {
    defined $def->{dataset_name} or die "Define dataset_name for bamplot first!";
    defined $def->{gene_names}   or die "Define gene_names for bamplot first, seperate by blank space!";
    defined $def->{add_chr}   or die "Define add_chr for bamplot first, check your genome sequence!";
  }

  my $config = {
    general => {
      task_name => $task,
      cluster   => $cluster
    },
    files => $def->{files}
  };

  if ( defined $def->{groups} ) {
    $config->{groups} = $def->{groups};
  }

  if ( defined $def->{pairs} ) {
    $config->{pairs} = $def->{pairs};
  }

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

  $config->{"fastqc_raw"} = {
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
  };
  $config->{"fastqc_raw_summary"} = {
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
  };
  push @individual, ("fastqc_raw");
  push @summary,    ("fastqc_raw_summary");

  my $configAlignment = {
    star => {
      class                     => "Alignment::STAR",
      perform                   => 1,
      target_dir                => $target_dir . "/star",
      option                    => "--twopassMode Basic",
      source_ref                => $source_ref,
      genome_dir                => $star_index,
      output_sort_by_coordinate => 1,
      sh_direct                 => 0,
      pbs                       => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "30gb"
      },
    },
    star_summary => {
      class      => "Alignment::STARSummary",
      perform    => 1,
      target_dir => $target_dir . "/star",
      option     => "",
      source_ref => [ "star", "_Log.final.out" ],
      sh_direct  => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    star_featurecount => {
      class      => "Count::FeatureCounts",
      perform    => 1,
      target_dir => $target_dir . "/star_featurecount",
      option     => "-g gene_id -t exon",
      source_ref => [ "star", "_Aligned.sortedByCoord.out.bam" ],
      gff_file   => $transcript_gtf,
      ispairend  => 1,
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    star_genetable => {
      class         => "CQS::CQSDatatable",
      perform       => 1,
      target_dir    => $target_dir . "/star_genetable",
      option        => "-k 0 -v 6 -e --fillMissingWithZero",
      source_ref    => "star_featurecount",
      name_map_file => $name_map_file,
      cqs_tools     => $cqstools,
      sh_direct     => 1,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    star_genetable_deseq2 => {
      class                => "Comparison::DESeq2",
      perform              => 1,
      target_dir           => $target_dir . "/star_genetable_deseq2",
      option               => "",
      source_ref           => "pairs",
      groups_ref           => "groups",
      countfile_ref        => "star_genetable",
      sh_direct            => 1,
      show_DE_gene_cluster => 0,
      pvalue               => 0.05,
      fold_change          => 2.0,
      min_median_read      => 5,
      use_pearson_in_hca   => $use_pearson_in_hca,
      pbs                  => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    star_genetable_correlation => {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      rCode                    => "usePearsonInHCA<-$use_pearson_in_hca; useGreenRedColorInHCA<-$use_green_red_color_in_hca; top25cvInHCA<-$top25cv_in_hca; ",
      target_dir               => $target_dir . "/star_genetable_correlation",
      rtemplate                => "countTableVisFunctions.R,countTableGroupCorrelation.R",
      output_file              => "parameterSampleFile1",
      output_file_ext          => ".Correlation.png",
      use_pearson_in_hca       => $def->{use_pearson_in_hca},
      parameterSampleFile1_ref => [ "star_genetable", ".count\$" ],
      parameterSampleFile2_ref => "groups",
      sh_direct                => 1,
      pbs                      => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "1",
        "mem"      => "10gb"
      },
    }
  };
  $config = merge( $config, $configAlignment );
  push @individual, ( "star", "star_featurecount" );
  push @summary, ( "star_summary", "star_genetable", "star_genetable_correlation", "star_genetable_deseq2" );

  if ( $def->{perform_rnaseqc} ) {
    $config->{rnaseqc} = {
      class          => "QC::RNASeQC",
      perform        => 1,
      target_dir     => "${target_dir}/rnaseqc",
      init_command   => $def->{rnaseqc_init_command},
      option         => "",
      source_ref     => [ "star", "_Aligned.sortedByCoord.out.bam" ],
      jar            => $def->{rnaseqc_jar},
      fasta_file     => $def->{fasta_file},
      rrna_fasta     => $def->{rrna_fasta},
      transcript_gtf => $transcript_gtf,
      pbs            => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    push( @summary, "rnaseqc" );
  }
  if ( $def->{perform_qc3bam} ) {
    $config->{star_qc3} = {
      class          => "QC::QC3bam",
      perform        => 1,
      target_dir     => $target_dir . "/star_qc3",
      option         => "",
      transcript_gtf => $transcript_gtf,
      qc3_perl       => $def->{qc3_perl},
      source_ref     => [ "star", "_Aligned.sortedByCoord.out.bam" ],
      pbs            => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    push( @summary, "star_qc3" );
  }

  if ( $def->{perform_bamplot} ) {
    $config->{gene_pos} = {
      class        => "Annotation::PrepareGenePosition",
      perform      => 1,
      target_dir   => $target_dir . "/gene_pos",
      option       => "",
      dataset_name => $def->{dataset_name},
      gene_names   => $def->{gene_names},
      add_chr      => $def->{add_chr},
      output_gff   => 1,
      pbs          => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      },
    };
    $config->{bamplot} = {
      class              => "Visualization::Bamplot",
      perform            => 1,
      target_dir         => "${target_dir}/bamplot",
      option             => "-g " . $def->{dataset_name} . " -y uniform -r --save-temp",
      source_ref         => [ "star", "_Aligned.sortedByCoord.out.bam" ],
      gff_file_ref       => "gene_pos",
      is_rainbow_color   => 0,
      is_single_pdf      => 0,
      is_draw_individual => 0,
      groups             => $def->{"plotgroups"},
      colors             => $def->{"colormaps"},
      sh_direct          => 1,
      pbs                => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "1",
        "mem"      => "10gb"
      },
    };
    push( @summary, "gene_pos", "bamplot" );
  }

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

sub performRNASeq {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getRNASeqConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

sub performRNASeqTask {
  my ( $def, $task ) = @_;

  my $config = getRNASeqConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
