#!/usr/bin/perl
package Pipeline::scRNASeq;

use strict;
use warnings;
use List::Util qw(first);
use File::Basename;
use Storable qw(dclone);
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(initializeScRNASeqDefaultOptions performScRNASeq performScRNASeqTask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeScRNASeqDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );

  initDefaultValue( $def, "perform_scRNABatchQC", 1 );

  initDefaultValue( $def, "perform_preprocessing", 0 );
  initDefaultValue( $def, "perform_mapping",       0 );
  initDefaultValue( $def, "perform_counting",      0 );
  initDefaultValue( $def, "perform_count_table",   0 );
  initDefaultValue( $def, "perform_correlation",   0 );
  initDefaultValue( $def, "perform_webgestalt",    0 );
  initDefaultValue( $def, "perform_report",        1 );
  initDefaultValue( $def, "perform_gsea",          0 );

  initDefaultValue( $def, "perform_cutadapt", 0 );

  initDefaultValue( $def, "featureCount_option",        "-g gene_id -t exon" );
  initDefaultValue( $def, "aligner",                    "star" );
  initDefaultValue( $def, "star_option",                "--twopassMode Basic --outSAMprimaryFlag AllBestScore" );
  initDefaultValue( $def, "use_pearson_in_hca",         1 );
  initDefaultValue( $def, "top25cv_in_hca",             0 );
  initDefaultValue( $def, "use_green_red_color_in_hca", 1 );
  initDefaultValue( $def, "output_bam_to_same_folder",  1 );
  initDefaultValue( $def, "show_label_PCA",             1 );

  initDefaultValue( $def, "max_thread",            8 );
  initDefaultValue( $def, "sequencetask_run_time", '24' );

  initDefaultValue( $def, "perform_keggprofile",      0 );
  initDefaultValue( $def, "keggprofile_useRawPValue", 0 );
  initDefaultValue( $def, "keggprofile_species",      "hsa" );
  initDefaultValue( $def, "keggprofile_pCut",         0.1 );

  initDefaultValue( $def, "is_paired_end",                   1 );
  initDefaultValue( $def, "DE_pvalue",                       0.05 );
  initDefaultValue( $def, "DE_use_raw_pvalue",               0 );
  initDefaultValue( $def, "DE_fold_change",                  2 );
  initDefaultValue( $def, "DE_export_significant_gene_name", 1 );
  initDefaultValue( $def, "DE_show_gene_cluster",            0 );
  initDefaultValue( $def, "DE_add_count_one",                0 );
  initDefaultValue( $def, "DE_top25only",                    0 );
  initDefaultValue( $def, "DE_detected_in_both_group",       0 );
  initDefaultValue( $def, "DE_perform_wilcox",               0 );
  initDefaultValue( $def, "DE_text_size",                    10 );
  initDefaultValue( $def, "DE_min_median_read",              5 );
  initDefaultValue( $def, "perform_DE_proteincoding_gene",   1 );
  initDefaultValue( $def, "perform_proteincoding_gene",      getValue( $def, "perform_DE_proteincoding_gene" ) );

  initDefaultValue( $def, "DE_outputPdf",  getValue( $def, "outputPdf",  0 ) );
  initDefaultValue( $def, "DE_outputPng",  getValue( $def, "outputPng",  1 ) );
  initDefaultValue( $def, "DE_outputTIFF", getValue( $def, "outputTIFF", 0 ) );
  initDefaultValue( $def, "DE_showVolcanoLegend", 1 );

  initDefaultValue( $def, "perform_seurat",      1 );
  initDefaultValue( $def, "Mtpattern",           "^MT-|^Mt-" );
  initDefaultValue( $def, "rRNApattern",         "^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]" );
  initDefaultValue( $def, "Remove_Mt_rRNA",      "FALSE" );
  initDefaultValue( $def, "nFeature_cutoff_min", 200 );
  initDefaultValue( $def, "nFeature_cutoff_max", 5000 );
  initDefaultValue( $def, "nCount_cutoff",       500 );
  initDefaultValue( $def, "mt_cutoff",           20 );
  initDefaultValue( $def, "resolution",          0.8 );
  initDefaultValue( $def, "pca_dims",            20 );
  initDefaultValue( $def, "details_rmd",         "" );

  initDefaultValue( $def, "by_integration", 0 );
  initDefaultValue( $def, "by_sctransform", 0 );

  initDefaultValue( $def, "perform_edgeR",      0 );
  initDefaultValue( $def, "DE_by_celltype",     1 );
  initDefaultValue( $def, "perform_webgestalt", 0 );

  return $def;
}

sub getScRNASeqConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeScRNASeqDefaultOptions($def);

  my $taskName = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $target_dir      = $def->{target_dir};
  my $groups_ref      = defined $def->{groups} ? "groups" : undef;
  my $aligner         = $def->{aligner};
  my $star_option     = $def->{star_option};
  my $count_table_ref = "files";

  if ( getValue( $def, "perform_scRNABatchQC" ) ) {
    $config->{scRNABatchQC} = {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      target_dir               => $target_dir . "/" . getNextFolderIndex($def) . "scRNABatchQC",
      rtemplate                => "../scRNA/scRNABatchQC.r",
      parameterSampleFile1_ref => "files",
      rCode                    => "webgestalt_organism='" . getValue( $def, "webgestalt_organism" ) . "'",
      sh_direct                => 1,
      pbs                      => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    };
    push( @$summary, "scRNABatchQC" );
  }

  if ( getValue( $def, "perform_seurat" ) ) {
    my $additional_rmd_files = "Functions.Rmd";

    my $seurat_name = "seurat" . ( getValue( $def, "by_sctransform" ) ? "_sct" : "" );
    $config->{$seurat_name} = {
      class                    => "CQS::UniqueRmd",
      perform                  => 1,
      target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $seurat_name,
      report_rmd_file          => "../scRNA/analysis.rmd",
      additional_rmd_files     => $additional_rmd_files,
      parameterSampleFile1_ref => "files",
      parameterSampleFile2     => {
        Mtpattern           => getValue( $def, "Mtpattern" ),
        rRNApattern         => getValue( $def, "rRNApattern" ),
        Remove_Mt_rRNA      => getValue( $def, "Remove_Mt_rRNA" ),
        nFeature_cutoff_min => getValue( $def, "nFeature_cutoff_min" ),
        nFeature_cutoff_max => getValue( $def, "nFeature_cutoff_max" ),
        nCount_cutoff       => getValue( $def, "nCount_cutoff" ),
        mt_cutoff           => getValue( $def, "mt_cutoff" ),
        resolution          => getValue( $def, "resolution" ),
        pca_dims            => getValue( $def, "pca_dims" ),
        species             => getValue( $def, "species" ),
        markers_file        => getValue( $def, "markers_file" ),
        details_rmd         => getValue( $def, "details_rmd" ),
        by_integration      => getValue( $def, "by_integration" ),
        by_sctransform      => getValue( $def, "by_sctransform" ),
        prefix              => $taskName,
      },
      output_file_ext => ".final.rds;.allmarkers.csv;.top10markers.csv;.cta_cluster.csv;_ur.html",
      sh_direct       => 1,
      pbs             => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    };
    push( @$summary, $seurat_name );

    my $countTaskname;
    if ( getValue( $def, "DE_by_celltype" ) ) {
      $countTaskname = $seurat_name . "_count_celltype";
    }
    else {
      $countTaskname = $seurat_name . "_count_cluster";
    }

    $config->{$countTaskname} = {
      class              => "CQS::UniqueR",
      perform            => 1,
      target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $countTaskname,
      rtemplate          => "../scRNA/scRNAcount.r",
      parameterFile1_ref => [ $seurat_name, ".final.rds" ],
      output_file_ext    => ".cluster.csv",
      rCode              => "DE_by_celltype=" . getValue( $def, "DE_by_celltype" ),
      sh_direct          => 1,
      pbs                => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    };
    push( @$summary, $countTaskname );

    if ( defined $def->{genes} ) {
      my $genesTaskname = $seurat_name . "_genes";
      $config->{$genesTaskname} = {
        class              => "CQS::UniqueR",
        perform            => 1,
        target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $genesTaskname,
        rtemplate          => "../scRNA/scRNAgenes.r",
        parameterFile1_ref => [ $seurat_name, ".final.rds" ],
        output_file_ext    => ".cluster.csv",
        rCode              => "genes='" . $def->{genes} . "'",
        sh_direct          => 1,
        pbs                => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary, $genesTaskname );
    }

    if ( getValue( $def, "perform_edgeR" ) ) {
      my $edgeRtaskname;
      if ( getValue( $def, "DE_by_celltype" ) ) {
        $edgeRtaskname = "edgeR_celltype";
      }
      else {
        $edgeRtaskname = "edgeR_cluster";
      }

      $config->{$edgeRtaskname} = {
        class                => "scRNA::scRNAedgeR",
        perform              => 1,
        target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $edgeRtaskname,
        rtemplate            => "../scRNA/scRNAedgeR.r",
        parameterFile1_ref   => [ $countTaskname, ".cluster.csv" ],
        parameterSampleFile1 => $def->{groups},
        parameterSampleFile2 => $def->{pairs},
        sh_direct            => 1,
        pbs                  => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary, $edgeRtaskname );

      if ( getValue( $def, "perform_webgestalt" ) ) {
        my $webgestaltTaskName = $edgeRtaskname . "_WebGestalt";
        $config->{$webgestaltTaskName} = {
          class            => "Annotation::WebGestaltR",
          perform          => 1,
          target_dir       => $target_dir . "/" . getNextFolderIndex($def) . $webgestaltTaskName,
          option           => "",
          source_ref       => [ $edgeRtaskname, "sig_genename.txt\$" ],
          organism         => getValue( $def, "webgestalt_organism" ),
          interestGeneType => $def->{interestGeneType},
          referenceSet     => $def->{referenceSet},
          sh_direct        => 1,
          pbs              => {
            "email"     => $email,
            "emailType" => $def->{emailType},
            "nodes"     => "1:ppn=1",
            "walltime"  => "23",
            "mem"       => "10gb"
          },
        };
        push @$summary, "$webgestaltTaskName";

        my $linkTaskName = $webgestaltTaskName . "_link_edgeR";
        $config->{$linkTaskName} = {
          class                      => "CQS::UniqueR",
          perform                    => 1,
          target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $linkTaskName,
          rtemplate                  => "../Annotation/WebGestaltDeseq2.r",
          rReportTemplate            => "../Annotation/WebGestaltDeseq2.rmd",
          output_to_result_directory => 1,
          output_perSample_file      => "parameterSampleFile1",
          output_perSample_file_ext  => ".html;.html.rds",
          parameterSampleFile1_ref   => [ $webgestaltTaskName, ".txt\$" ],
          parameterSampleFile2_ref   => [ $edgeRtaskname, "sig.csv\$" ],
          sh_direct                  => 1,
          rCode                      => "",
          pbs                        => {
            "email"     => $def->{email},
            "emailType" => $def->{emailType},
            "nodes"     => "1:ppn=1",
            "walltime"  => "23",
            "mem"       => "10gb"
          },
        };
        push( @$summary, $linkTaskName );
      }
    }
  }

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => $individual,
      step2 => $summary,
    },
    sh_direct => 0,
    cluster   => $cluster,
    pbs       => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => $def->{sequencetask_run_time},
      "mem"       => "40gb"
    },
  };

  return ($config);
}

sub performScRNASeq {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getScRNASeqConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

sub performScRNASeqTask {
  my ( $def, $task ) = @_;

  my $config = performScRNASeq($def);

  performTask( $config, $task );

  return $config;
}

1;
