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
use Storable qw(dclone);
use scRNA::Modules;

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

  initDefaultValue( $def, "perform_enclone_only",      0 );

  initDefaultValue( $def, "perform_seurat",      1 );
  initDefaultValue( $def, "Mtpattern",           "^MT-|^Mt-" );
  initDefaultValue( $def, "rRNApattern",         "^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]" );
  initDefaultValue( $def, "Remove_Mt_rRNA",      "FALSE" );
  initDefaultValue( $def, "nFeature_cutoff_min", 200 );
  initDefaultValue( $def, "nFeature_cutoff_max", 10000 );
  initDefaultValue( $def, "nCount_cutoff",       500 );
  initDefaultValue( $def, "mt_cutoff",           20 );
  initDefaultValue( $def, "resolution",          0.5 );
  initDefaultValue( $def, "pca_dims",            20 );
  initDefaultValue( $def, "details_rmd",         "" );

  initDefaultValue( $def, "by_integration",        0 );
  initDefaultValue( $def, "by_sctransform",        0 );
  initDefaultValue( $def, "pool_sample",           0 );
  initDefaultValue( $def, "batch_for_integration", 0 );

  initDefaultValue( $def, "perform_recluster",      0 );
  initDefaultValue( $def, "perform_rename_cluster", 0 );
  initDefaultValue( $def, "perform_antibody_vis",   0 );
  initDefaultValue( $def, "antibody_pattern",       "^CD" );

  initDefaultValue( $def, "perform_edgeR", 0 );

  initDefaultValue( $def, "DE_by_celltype", 1 );
  initDefaultValue( $def, "DE_by_cluster",  0 );
  if ( ( not getValue( $def, "DE_by_celltype" ) ) && ( not getValue( $def, "DE_by_cluster" ) ) ) {
    initDefaultValue( $def, "DE_by_celltype", 1 );
  }

  initDefaultValue( $def, "DE_by_sample", 0 );
  initDefaultValue( $def, "DE_by_cell",   1 );
  if ( ( not getValue( $def, "DE_by_sample" ) ) && ( not getValue( $def, "DE_by_sample" ) ) ) {
    initDefaultValue( $def, "DE_by_cell", 1 );
  }

  initDefaultValue( $def, "DE_by_cell_filter_minTPM",         1 );
  initDefaultValue( $def, "DE_by_cell_filter_cellPercentage", 0.25 );

  initDefaultValue( $def, "DE_by_sample_filter_minTPM",         1 );
  initDefaultValue( $def, "DE_by_sample_filter_cellPercentage", 0.5 );

  initDefaultValue( $def, "perform_webgestalt", 0 );

  initDefaultValue( $def, "perform_CHETAH", 0 );
  
  return $def;
}

sub addAntibodyTask {
  my ( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name ) = @_;

  my $pattern = getValue( $def, "antibody_pattern" );

  my $taskname = $cluster_task_name . "_antibody_vis";
  $config->{$taskname} = {
    class              => "CQS::UniqueR",
    perform            => 1,
    target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $taskname,
    rtemplate          => "../scRNA/scRNAantibody.r",
    parameterFile1_ref => [ $seurat_name, ".final.rds" ],
    parameterFile3_ref => [ $cluster_task_name, $cluster_file ],
    output_file_ext    => ".cluster.csv",
    rCode              => "celltype_name='$celltype_name'; cluster_name='$cluster_name'; antibody_pattern='$pattern'; ",
    sh_direct          => 1,
    pbs                => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push( @$summary, $taskname );
}

sub addMarkerGenes {
  my ( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name, $marker_name, $marker_file, $samples ) = @_;
  
  my $markerGenesTaskname = $cluster_task_name . "_" . $marker_name;
  $config->{$markerGenesTaskname} = {
    class              => "CQS::UniqueR",
    perform            => 1,
    target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $markerGenesTaskname,
    rtemplate          => "../scRNA/scRNAMarkerGenes.r",
    parameterFile1_ref => [ $seurat_name, ".final.rds" ],
    parameterFile2     => $marker_file,
    parameterFile3_ref => [ $cluster_task_name, $cluster_file ],
    output_file_ext    => ".cluster.csv",
    rCode              => "celltype_name='$celltype_name'; cluster_name='$cluster_name';",
    sh_direct          => 1,
    pbs                => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  
  if (defined $samples){
    $config->{$markerGenesTaskname}{rCode} = $config->{$markerGenesTaskname}{rCode} . "samples='$samples';";
  }
  push( @$summary, $markerGenesTaskname );
}
sub addGeneTask {
  my ( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name ) = @_;

  my $marker_genes = $def->{plot_marker_genes};
  if (not defined $marker_genes){
    $marker_genes = {};
  }
  
  if ( defined $def->{marker_genes_file} ) {
    $marker_genes->{"marker"} = {
      file => $def->{marker_genes_file}
    };
  }
  
  if (defined $def->{pathway_genes_file}){
    $marker_genes->{"pathway"} = {
      file => $def->{pathway_genes_file}
    };
  }

  for my $key (sort keys %$marker_genes){
    my $file = $marker_genes->{$key}{file};
    my $samples = $marker_genes->{$key}{samples};
    addMarkerGenes($config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name, "genes_" . $key, $file, $samples);
  }


  if ( defined $def->{genes} ) {
    my $dotPlotOnly = getValue($def, "genesDotPlotOnly", "0");
    my $genes = $def->{genes};
    $genes =~ s/\n/;/g;
    $genes =~ s/\s/;/g;
    $genes =~ s/;;/;/g;
    my $genesTaskname = $cluster_task_name . "_genes";
    $config->{$genesTaskname} = {
      class              => "CQS::UniqueR",
      perform            => 1,
      target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $genesTaskname,
      rtemplate          => "../scRNA/scRNAgenes.r",
      parameterFile1_ref => [ $seurat_name, ".final.rds" ],
      parameterFile3_ref => [ $cluster_task_name, $cluster_file ],
      output_file_ext    => ".cluster.csv",
      rCode              => "genes='" . $genes . "'; celltype_name='$celltype_name'; cluster_name='$cluster_name'; dotPlotOnly=$dotPlotOnly;",
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
}

sub addDeseq2BySampleTask {
  my ( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name, $bBetweenCluster, $DE_by_celltype) = @_;
  my $rCode = "pvalue=" . getValue( $def, "DE_pvalue" ) . ";useRawPvalue=" . getValue( $def, "DE_use_raw_pvalue" ) . ";foldChange=" . getValue( $def, "DE_fold_change" );

  #die "Found";

  my $prefix  = $cluster_task_name;
  my $curClusterName = undef;
  my $curClusterDisplayName = undef;
  if ($DE_by_celltype) {
    $curClusterName = $celltype_name;
    $prefix  = $prefix . "_inCelltype";
  }
  else {
    $curClusterName = $cluster_name;
    $prefix  = $prefix . "_inCluster";
  }

  $rCode         = $rCode . ";DE_by_cell=0;filter_minTPM=" . getValue( $def, "DE_by_sample_filter_minTPM" ) . ";filter_samplePercentage=" . getValue( $def, "DE_by_sample_filter_cellPercentage" );
  $prefix = $prefix . "_bySample";

  $rCode  = $rCode . ";bBetweenCluster=0";
  my $groups = getValue( $def, "groups" );
  my $pairs  = getValue( $def, "pairs" );
  $rCode = $rCode . ";cluster_name='" . $curClusterName . "'";

  my $deseq2table_taskname=$prefix . "_table";
  $config->{$deseq2table_taskname} = {
    class                => "CQS::UniqueR",
    perform              => 1,
    target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $deseq2table_taskname,
    rtemplate            => "../scRNA/deseq2table.r",
    parameterFile1_ref   => [ $seurat_name, ".final.rds" ],
    parameterFile2_ref   => [ $cluster_task_name, $cluster_file ],
    parameterSampleFile1 => $groups,
    parameterSampleFile2 => $pairs,
    output_file_ext      => ".deseq2.define.txt",
    rCode                => $rCode,
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push( @$summary, $deseq2table_taskname );
}

sub addEdgeRTask {
  my ( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name, $bBetweenCluster, $DE_by_celltype, $DE_by_cell ) = @_;
  my $rCode = "pvalue=" . getValue( $def, "DE_pvalue" ) . ";useRawPvalue=" . getValue( $def, "DE_use_raw_pvalue" ) . ";foldChange=" . getValue( $def, "DE_fold_change" );

  my $edgeRtaskname  = $cluster_task_name . "_edgeR";
  my $groups         = undef;
  my $pairs          = undef;
  my $curClusterName = undef;
  my $curClusterDisplayName = undef;
  if ($bBetweenCluster) {
    $edgeRtaskname  = $edgeRtaskname . "_betweenCluster_byCell";
    $curClusterName = getValue( $def, "DE_cluster_name" );
    $curClusterDisplayName = getValue( $def, "DE_cluster_display_name", $curClusterName );
    $rCode  = $rCode . ";filter_minTPM=" . getValue( $def, "DE_by_cell_filter_minTPM" ) . ";filter_cellPercentage=" . getValue( $def, "DE_by_cell_filter_cellPercentage" ) . ";bBetweenCluster=1";
    $groups = getValue( $def, "DE_cluster_groups" );
    $pairs  = getValue( $def, "DE_cluster_pairs" );
  }
  else {
    if ($DE_by_celltype) {
      $curClusterName = $celltype_name;
      $edgeRtaskname  = $edgeRtaskname . "_inCelltype";
    }
    else {
      $curClusterName = $cluster_name;
      $edgeRtaskname  = $edgeRtaskname . "_inCluster";
    }

    if ($DE_by_cell) {
      $rCode         = $rCode . ";DE_by_cell=1;filter_minTPM=" . getValue( $def, "DE_by_cell_filter_minTPM" ) . ";filter_cellPercentage=" . getValue( $def, "DE_by_cell_filter_cellPercentage" );
      $edgeRtaskname = $edgeRtaskname . "_byCell";
    }
    else {
      $rCode         = $rCode . ";DE_by_cell=0;filter_minTPM=" . getValue( $def, "DE_by_sample_filter_minTPM" ) . ";filter_samplePercentage=" . getValue( $def, "DE_by_sample_filter_cellPercentage" );
      $edgeRtaskname = $edgeRtaskname . "_bySample";
    }

    $rCode  = $rCode . ";bBetweenCluster=0";
    $groups = getValue( $def, "groups" );
    $pairs  = getValue( $def, "pairs" );
  }
  $rCode = $rCode . ";cluster_name='" . $curClusterName . "'";

  $config->{$edgeRtaskname} = {
    class                => "CQS::UniqueR",
    perform              => 1,
    target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $edgeRtaskname,
    rtemplate            => "../scRNA/edgeR.r",
    parameterFile1_ref   => [ $seurat_name, ".final.rds" ],
    parameterFile2_ref   => [ $cluster_task_name, $cluster_file ],
    parameterSampleFile1 => $groups,
    parameterSampleFile2 => $pairs,
    output_file_ext      => ".edgeR.files.csv",
    rCode                => $rCode,
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push( @$summary, $edgeRtaskname );

  my $vistaskname = $edgeRtaskname . "_vis";
  $config->{$vistaskname} = {
    class              => "CQS::UniqueR",
    perform            => 1,
    target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $vistaskname,
    rtemplate          => "../scRNA/edgeRvis.r",
    parameterFile1_ref => [ $seurat_name, ".final.rds" ],
    parameterFile2_ref => [$edgeRtaskname],
    parameterFile3_ref => [ $cluster_task_name, $cluster_file ],
    output_file_ext    => ".edgeRvis.files.csv",
    rCode              => "cluster_name='" . $curClusterName . "';bBetweenCluster=" . $bBetweenCluster,
    sh_direct          => 1,
    pbs                => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push( @$summary, $vistaskname );

  if($bBetweenCluster) {
    my $vistaskname2 = $edgeRtaskname . "_dotplot";
    $config->{$vistaskname2} = {
      class              => "CQS::UniqueR",
      perform            => 1,
      target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $vistaskname2,
      rtemplate          => "../scRNA/scRNA_func.r;../scRNA/edgeRdotplot.r",
      parameterFile1_ref => [ $seurat_name, ".final.rds" ],
      parameterFile2_ref => [$edgeRtaskname],
      parameterFile3_ref => [ $cluster_task_name, $cluster_file ],
      parameterSampleFile1 => {
        cluster_name => getValue( $def, "DE_clusters_name", $curClusterName ),
        display_cluster_name => getValue( $def, "DE_clusters_display_name", $curClusterDisplayName ),
        gene_number => getValue( $def, "DE_dotplot_gene_number", 20 ),
      },
      output_file_ext    => ".edgeRvis2.files.csv",
      rCode              => "",
      sh_direct          => 1,
      pbs                => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    };
    push( @$summary, $vistaskname2 );
  }

  if ( getValue( $def, "perform_webgestalt" ) ) {
    my $webgestaltTaskName = $edgeRtaskname . "_WebGestalt";
    $config->{$webgestaltTaskName} = {
      class              => "CQS::UniqueR",
      perform            => 1,
      target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $webgestaltTaskName,
      rtemplate          => "../Annotation/WebGestaltReportFunctions.r;../Annotation/WebGestaltR_all.r",
      parameterFile1_ref => [$edgeRtaskname],
      parameterSampleFile1 => {
        organism         => getValue( $def, "webgestalt_organism" ),
        interestGeneType => $def->{interestGeneType},
        referenceSet     => $def->{referenceSet},
      },
      output_file_ext    => ".WebGestaltR.files.csv",
      rCode              => "",
      sh_direct          => 1,
      pbs                => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    };
    push @$summary, "$webgestaltTaskName";

    my $linkTaskName = $webgestaltTaskName . "_link_edgeR";
    $config->{$linkTaskName} = {
      class                      => "CQS::UniqueR",
      perform                    => 1,
      target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $linkTaskName,
      rtemplate                  => "../Annotation/WebGestaltReportFunctions.r;../Annotation/WebGestaltDeseq2_all.r",
      rReportTemplate            => "../Annotation/WebGestaltDeseq2.rmd",
      output_to_result_directory => 1,
      parameterFile1_ref   => [ $webgestaltTaskName ],
      parameterFile2_ref   => [ $edgeRtaskname ],
      output_file_ext    => ".link.files.csv",
      sh_direct                  => 1,
      rCode                      => "",
      pbs                        => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "23",
        "mem"       => "10gb"
      },
    };
    push( @$summary, $linkTaskName );
  }

  if ( getValue( $def, "perform_gsea" ) ) {
    my $gsea_jar        = $def->{gsea_jar}        or die "Define gsea_jar at definition first";
    my $gsea_db         = $def->{gsea_db}         or die "Define gsea_db at definition first";
    my $gsea_categories = $def->{gsea_categories} or die "Define gsea_categories at definition first";
    my $gsea_makeReport = 1;
    if ( defined $def->{gsea_makeReport} ) {
      $gsea_makeReport = $def->{gsea_makeReport};
    }

    my $gseaTaskName = $edgeRtaskname . "_GSEA";

    #my $gseaCategories = "'h.all.v6.1.symbols.gmt','c2.all.v6.1.symbols.gmt','c5.all.v6.1.symbols.gmt','c6.all.v6.1.symbols.gmt','c7.all.v6.1.symbols.gmt'";
    $config->{$gseaTaskName} = {
      class                      => "CQS::UniqueR",
      perform                    => 1,
      target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $gseaTaskName,
      rtemplate                  => "GSEAPerform.R",
      rReportTemplate            => "GSEAReport.Rmd",
      output_to_result_directory => 1,
      output_file_ext            => ".gsea.files.csv",
      parameterFile1_ref         => [ $edgeRtaskname, ".edgeR.files.csv\$" ],
      sh_direct                  => 1,
      rCode                      => "gseaDb='" . $gsea_db . "'; gseaJar='" . $gsea_jar . "'; gseaCategories=c(" . $gsea_categories . "); makeReport=" . $gsea_makeReport . ";",
      pbs                        => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "23",
        "mem"       => "10gb"
      },
    };
    push( @$summary, $gseaTaskName );

  #   if($bBetweenCluster) {
  #     my @gsea_report_files = ();
  #     my @gsea_report_names = ();
  #     my $pairs = $config->{pairs};
  #     for my $key ( keys %$pairs ) {
  #       push( @gsea_report_files, $gseaTaskName, "/.*." . $key . ".edgeR_GSEA.rnk.gsea.csv" );
  #       push( @gsea_report_names, "gsea_" . $key );
  #     }

  #     my $gsea_report = $gseaTaskName . "_report";
  #     $config->{$gsea_report} = {
  #       class                      => "CQS::BuildReport",
  #       perform                    => 1,
  #       target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $gsea_report,
  #       report_rmd_file            => "GSEAReport.Rmd",
  #       additional_rmd_files       => "../Pipeline/Pipeline.Rmd;Functions.Rmd",
  #       parameterSampleFile1_ref   => \@gsea_report_files,
  #       parameterSampleFile1_names => \@gsea_report_names,
  #       parameterSampleFile3       => [],
  #       sh_direct                  => 1,
  #       pbs                        => {
  #         "nodes"     => "1:ppn=1",
  #         "walltime"  => "1",
  #         "mem"       => "10gb"
  #       },
  #     };
  #     push( @$summary, $gsea_report );
  #   }
  }
  
  return ($edgeRtaskname);
}

sub getScRNASeqConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeScRNASeqDefaultOptions($def);

  my $project_name = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $target_dir      = $def->{target_dir};
  my $groups_ref      = defined $def->{groups} ? "groups" : undef;
  my $aligner         = $def->{aligner};
  my $star_option     = $def->{star_option};
  my $count_table_ref = "files";

  $config->{bam_files} = $def->{bam_files};

  my $perform_split_hto_samples = getValue($def, "perform_split_hto_samples", 0);
  my $has_vdj_json_files = defined $def->{vdj_json_files};

  my $clonotype_4_convert;
  if (defined $def->{vdj_json_files}){
    if ((not defined $def->{files}) || (not $perform_split_hto_samples)) {
    $config->{vdj_json_files} = $def->{vdj_json_files};
      addClonotypeMerge($config, $def, $summary, $target_dir, "clonotype_1_merge", ["vdj_json_files", "all_contig_annotations.json"]);
      addEnclone($config, $def, $summary, "clonotype_2_enclone", $target_dir, ["clonotype_1_merge", ".json\$"] );
      $clonotype_4_convert = addEncloneToClonotype($config, $def, $summary, $target_dir, "clonotype_3_convert", "clonotype_2_enclone", ["clonotype_1_merge", ".cdr3\$"]);
    }
  }

  if (defined $def->{files}){
    my @report_files = ();
    my @report_names = ();
    my $hto_name = undef;
    my $hto_ref = undef;
    my $hto_sample_file = undef;
    my $hto_summary = undef;
    if( $perform_split_hto_samples ) {
      my $r_script = undef;
      my $folder = undef;
      if ( getValue($def, "split_hto_samples_by_cutoff", 0) ) {
        $r_script = "../scRNA/split_samples_cutoff.r";
        $folder = "hto_samples_cutoff";
      } else {
        $r_script = "../scRNA/split_samples.r";
        $folder = "hto_samples_HTODemux";
      }

      my $files = $def->{files};
      my $hto_file_names = $def->{hto_file_names};
      my $hto_file_ref = "files";
      if(defined $hto_file_names){
        my $hto_files = {};
        for my $hto_name (@$hto_file_names){
          $hto_files->{$hto_name} = $files->{$hto_name};
        }
        $config->{hto_files} = $hto_files;
        $hto_file_ref = "hto_files";
        #print("hto_files=" . Dumper($hto_files));
      }
      #print("hto_file_ref=" . $hto_file_ref . "\n");
      $hto_name = "hto_samples";
      $config->{$hto_name} = {
        class => "CQS::ProgramWrapperOneToOne",
        target_dir => "${target_dir}/$folder",
        interpretor => getValue($def, "R", "R") . " --vanilla -f ",
        program => $r_script,
        check_program => 1,
        option => "--args __FILE__ __OUTPUT__ " . getValue($def, "hto_regex", ""),
        source_arg => "",
        source_ref => $hto_file_ref,
        output_arg => "",
        output_file_prefix => ".HTO",
        output_file_ext => ".HTO.class.dist.png;.HTO.csv",
        output_to_same_folder => 0,
        can_result_be_empty_file => 0,
        sh_direct   => 1,
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$individual, "hto_samples" );

      $hto_ref = [ $hto_name, ".HTO.csv" ];

      $hto_summary = $hto_name . "_summary";
      $config->{$hto_summary} = {
        class => "CQS::UniqueR",
        target_dir => "${target_dir}/${folder}_summary",
        rtemplate => "../scRNA/split_samples_summary.r",
        option => "",
        parameterSampleFile1_ref => $hto_ref,
        parameterSampleFile2 => $def->{"HTO_name_map"},
        output_file => "",
        output_file_ext => ".HTO.summary.csv;.HTO.summary.global.png",
        sh_direct   => 1,
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      },

      push( @$individual, $hto_summary );

      push (@report_files, ($hto_summary, ".HTO.summary.global.png"));
      push (@report_names, "hto_summary_png");

      if(defined $def->{HTO_samples}){
        $hto_sample_file = write_HTO_sample_file($def);
      }

      if(defined $def->{bam_files}){
        if (not defined $def->{HTO_samples}) {
          die "Define HTO_samples for split bam files";
        }

        $config->{HTO_samples} = $def->{HTO_samples};
        $config->{bam_files} = $def->{bam_files};
        $config->{"hto_bam"} = {
          class => "CQS::ProgramWrapperOneToOne",
          target_dir => "${target_dir}/hto_bam",
          interpretor => "python3",
          program => "../scRNA/split_samples.py",
          check_program => 1,
          option => "-o .",
          source_arg => "-i",
          source_ref => $hto_ref,
          parameterSampleFile2_arg => "-b",
          parameterSampleFile2_ref => ["bam_files"],
          parameterSampleFile3_arg => "-s",
          parameterSampleFile3_ref => ["HTO_samples"],
          output_arg => "-o",
          output_file_prefix => "",
          output_file_ext => ".bam",
          output_to_same_folder => 1,
          can_result_be_empty_file => 0,
          sh_direct   => 1,
          pbs => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        };
        push( @$individual, "hto_bam" );

        addArcasHLA($config, $def, $individual, $target_dir, $project_name, "hto_bam", "hto_bam");        
      }

      if(defined $def->{vdj_json_files}){
        if (not defined $def->{HTO_samples}) {
          die "Define HTO_samples for split vdj json files";
        }

        $config->{HTO_samples} = $def->{HTO_samples};
        $config->{vdj_json_files} = $def->{vdj_json_files};
        $config->{"hto_clonotype_1_split"} = {
          class => "CQS::ProgramWrapperOneToManyFile",
          target_dir => "${target_dir}/hto_clonotype_1_split",
          interpretor => "python3",
          program => "../scRNA/clonotype_split.py",
          check_program => 1,
          option => "-o .",
          source_arg => "-i",
          source_ref => "vdj_json_files",
          parameterSampleFile2_arg => "-c",
          parameterSampleFile2_ref => $hto_ref,
          parameterSampleFile3_arg => "-s",
          parameterSampleFile3_ref => ["HTO_samples"],
          output_file => "parameterSampleFile3",
          output_arg => "-o",
          output_file_prefix => "",
          output_file_ext => "all_contig_annotations.json",
          output_to_same_folder => 0,
          samplename_in_result => 0,
          output_file_key => 0,
          can_result_be_empty_file => 0,
          sh_direct   => 1,
          pbs => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        };
        push( @$individual, "hto_clonotype_1_split" );

        addClonotypeMerge($config, $def, $summary, $target_dir, "hto_clonotype_2_merge", ["hto_clonotype_1_split", "all_contig_annotations.json"]);
        addEnclone($config, $def, $summary, "hto_clonotype_3_enclone", $target_dir, ["hto_clonotype_2_merge", ".json\$"] );
        $clonotype_4_convert = addEncloneToClonotype($config, $def, $summary, $target_dir, "hto_clonotype_4_convert", "hto_clonotype_3_enclone", ["hto_clonotype_2_merge", ".cdr3\$"]);
      }
    }else{
      if(defined $def->{bam_files}){
        addArcasHLA($config, $def, $individual, $target_dir, $project_name, "", "bam_files");        
      }
    }

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

    my $seurat_name;
    if ( getValue( $def, "perform_seurat" ) ) {
      $seurat_name = "seurat" . ( getValue( $def, "by_sctransform" ) ? "_sct" : "" ) . ( getValue( $def, "by_integration" ) ? "_igr" : "" ) . ( getValue( $def, "pool_sample" ) ? "_pool" : "" );
      $config->{$seurat_name} = {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $seurat_name,
        rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_cluster.r",
        parameterSampleFile1_ref => "files",
        parameterSampleFile2     => {
          Mtpattern             => getValue( $def, "Mtpattern" ),
          rRNApattern           => getValue( $def, "rRNApattern" ),
          Remove_Mt_rRNA        => getValue( $def, "Remove_Mt_rRNA" ),
          nFeature_cutoff_min   => getValue( $def, "nFeature_cutoff_min" ),
          nFeature_cutoff_max   => getValue( $def, "nFeature_cutoff_max" ),
          nCount_cutoff         => getValue( $def, "nCount_cutoff" ),
          mt_cutoff             => getValue( $def, "mt_cutoff" ),
          species               => getValue( $def, "species" ),
          resolution            => getValue( $def, "resolution" ),
          pca_dims              => getValue( $def, "pca_dims" ),
          by_integration        => getValue( $def, "by_integration" ),
          by_sctransform        => getValue( $def, "by_sctransform" ),
          pool_sample           => getValue( $def, "pool_sample" ),
          batch_for_integration => getValue( $def, "batch_for_integration" ),
          hto_sample_file       => $hto_sample_file,
        },
        parameterSampleFile3 => $def->{"batch_for_integration_groups"},
        parameterSampleFile4 => $def->{"pool_sample_groups"},
        parameterSampleFile5_ref => $hto_ref,
        output_file_ext      => ".final.rds",
        output_other_ext  => ".cluster.csv;.allmarkers.csv;.top10markers.csv;.cluster.normByUpQuantile.csv",
        sh_direct            => 1,
        pbs                  => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary, $seurat_name );

      my $celltype=$seurat_name . "_celltype";

      if(getValue( $def, "annotate_tcell", 0)){
        getValue( $def, "HLA_panglao5_file");
        getValue( $def, "tcell_markers_file");
      }

      $config->{$celltype} = {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $celltype,
        rtemplate                => "../scRNA/scRNA_func.r,../scRNA/celltype_annotation.r",
        parameterFile1_ref =>  [$seurat_name, ".cluster.normByUpQuantile.csv"],
        parameterFile2_ref =>  [$seurat_name, ".cluster.csv"],
        parameterSampleFile1    => {
          species               => getValue( $def, "species" ),
          db_markers_file       => getValue( $def, "markers_file" ),
          curated_markers_file  => getValue( $def, "curated_markers_file", "" ),
          annotate_tcell        => getValue( $def, "annotate_tcell", 0),
          remove_subtype        => getValue( $def, "remove_subtype", ""),
          HLA_panglao5_file     => getValue( $def, "HLA_panglao5_file", "" ),
          tcell_markers_file    => getValue( $def, "tcell_markers_file", ""),
        },
        output_file_ext      => ".celltype.csv",
        output_other_ext  => ".celltype_cluster.csv;.celltype.rds",
        sh_direct            => 1,
        pbs                  => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary, $celltype );
      my $cluster_task_name = $celltype;
      my $celltype_file     = ".celltype.csv";
      my $cluster_file      = ".celltype_cluster.csv";
      my $celltype_name     = "cellactivity_clusters";
      my $cluster_name      = "seurat_cellactivity_clusters";

      push(@report_files, ($seurat_name, ".final.rds", $celltype, ".celltype.csv", $celltype, ".celltype.rds"));
      push(@report_names, ("seurat_rds", "activity_celltype", "activity_rds"));

      # my $additional_rmd_files = "Functions.Rmd";
      # $config->{$seurat_name} = {
      #   class                    => "CQS::UniqueRmd",
      #   perform                  => 1,
      #   target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $seurat_name,
      #   report_rmd_file          => "../scRNA/analysis.rmd",
      #   additional_rmd_files     => $additional_rmd_files,
      #   parameterSampleFile1_ref => "files",
      #   parameterSampleFile2     => {
      #     Mtpattern             => getValue( $def, "Mtpattern" ),
      #     rRNApattern           => getValue( $def, "rRNApattern" ),
      #     Remove_Mt_rRNA        => getValue( $def, "Remove_Mt_rRNA" ),
      #     nFeature_cutoff_min   => getValue( $def, "nFeature_cutoff_min" ),
      #     nFeature_cutoff_max   => getValue( $def, "nFeature_cutoff_max" ),
      #     nCount_cutoff         => getValue( $def, "nCount_cutoff" ),
      #     mt_cutoff             => getValue( $def, "mt_cutoff" ),
      #     resolution            => getValue( $def, "resolution" ),
      #     pca_dims              => getValue( $def, "pca_dims" ),
      #     species               => getValue( $def, "species" ),
      #     markers_file          => getValue( $def, "markers_file" ),
      #     annotate_tcell        => getValue( $def, "annotate_tcell", 0),
      #     HLA_panglao5_file     => getValue( $def, "HLA_panglao5_file", "" ),
      #     tcell_markers_file    => getValue( $def, "tcell_markers_file", ""),
      #     details_rmd           => getValue( $def, "details_rmd" ),
      #     by_integration        => getValue( $def, "by_integration" ),
      #     by_sctransform        => getValue( $def, "by_sctransform" ),
      #     pool_sample           => getValue( $def, "pool_sample" ),
      #     batch_for_integration => getValue( $def, "batch_for_integration" ),
      #     hto_sample_file       => $hto_sample_file,
      #     prefix                => $project_name,
      #   },
      #   parameterSampleFile3 => $def->{"batch_for_integration_groups"},
      #   parameterSampleFile4 => $def->{"pool_sample_groups"},
      #   parameterSampleFile5_ref => $hto_ref,
      #   output_file_ext      => ".final.rds",
      #   output_other_ext  => ".cluster.csv;.allmarkers.csv;.top10markers.csv;.cluster.normByUpQuantile.csv;_ur.html",
      #   sh_direct            => 1,
      #   pbs                  => {
      #     "email"     => $def->{email},
      #     "emailType" => $def->{emailType},
      #     "nodes"     => "1:ppn=1",
      #     "walltime"  => "1",
      #     "mem"       => "10gb"
      #   },
      # };
      # push( @$summary, $seurat_name );

      my $find_markers = $seurat_name . "_celltype_markers";
      $config->{$find_markers} = {
        class                => "CQS::UniqueR",
        perform              => 1,
        target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $find_markers,
        rtemplate            => "../scRNA/scRNA_func.r,../scRNA/celltype_markers.r",
        parameterFile1_ref   => [ $seurat_name, ".final.rds" ],
        parameterFile2_ref   => [ $seurat_name, ".cluster.csv" ],
        parameterFile3_ref   => [ $celltype, ".celltype.csv" ],
        parameterSampleFile1 => {
          by_sctransform        => getValue( $def, "by_sctransform" ),
          celltype_name         => $celltype_name
        },
        output_file_ext      => ".all_display_markers.csv",
        output_other_ext     => ".all_in_markers.csv",
        sh_direct            => 1,
        pbs                  => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary, $find_markers );

      if (getValue( $def, "perform_scMRMA", 0 ) ) {
        addScMRMA( $config, $def, $summary, $target_dir, $project_name, $seurat_name );
      }

      my $parameterSampleFile5_ref = undef;
      if (getValue( $def, "perform_CHETAH", 0 ) ) {
        my $CHETAH_name= $seurat_name . "_CHETAH";
        addCHETAH( $config, $def, $summary, $target_dir, $project_name, $CHETAH_name, $seurat_name );
        push @report_files, ($CHETAH_name, ".CHETAH.png");
        push @report_names, "chetah_png";
      }

      my $parameterSampleFile6_ref = undef;
      if (getValue( $def, "perform_signac", 0 ) ) {
        my $signac_name = $seurat_name . "_signac";
        addSignac( $config, $def, $summary, $target_dir, $project_name, $signac_name, $seurat_name );
        push @report_files, ($signac_name, ".signac.png");
        push @report_names, "signac_png";
      }

      addGeneTask( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name );

      if ( getValue( $def, "perform_rename_cluster" ) ) {
        my $rename_cluster = $celltype . "_rename";
        $config->{$rename_cluster} = {
          class                => "CQS::UniqueR",
          perform              => 1,
          target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $rename_cluster,
          rtemplate            => "../scRNA/scRNA_func.r,../scRNA/rename_cluster.r",
          parameterFile1_ref   => [ $seurat_name, ".final.rds" ],
          parameterFile2_ref   => [ $celltype, ".cluster.csv" ],
          parameterFile3_ref   => [ $celltype, ".celltype.csv" ],
          parameterSampleFile1 => {
            celltype_name         => $celltype_name
          },
          parameterSampleFile2 => $def->{rename_cluster},
          output_file_ext      => ".rename_celltype.csv",
          output_other_ext     => ".rename_cluster.csv;.rename_cluster.png;.rename_cluster.summery.csv",
          sh_direct            => 1,
          pbs                  => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        };
        push( @$summary, $rename_cluster );

        push @report_files, ($rename_cluster, ".rename_cluster.png");
        push @report_names, "rename_png";

        $cluster_task_name = $rename_cluster;
        $celltype_file     = ".rename_celltype.csv";
        $cluster_file      = ".rename_cluster.csv";
        $celltype_name     = "renamed_cellactivity_clusters";
        $cluster_name      = "seurat_renamed_cellactivity_clusters";

        addGeneTask( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name );
      }

      my $report_name= "report";
      my $additional_rmd_files = "Functions.Rmd;../scRNA/scRNA_func.r";
      $config->{$report_name} = {
        class                    => "CQS::BuildReport",
        perform                  => 1,
        target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $report_name,
        report_rmd_file          => "../scRNA/report.rmd",
        additional_rmd_files     => $additional_rmd_files,
        parameterSampleFile1_ref => \@report_files,
        parameterSampleFile1_names => \@report_names,
        parameterSampleFile2 => merge({
            prefix => $project_name,
            summary_layer_file => $def->{summary_layer_file},
            celltype_name => $celltype_name
          }, merge($config->{$seurat_name}{parameterSampleFile2}, $config->{$celltype}{parameterSampleFile1})),
        parameterSampleFile3 => [],
        output_file_ext      => "_ur.html",
        sh_direct            => 1,
        pbs                  => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      if(defined $hto_name){
        $config->{$report_name}{parameterSampleFile4_ref} = [$hto_name, ".HTO.class.dist.png"];
      }
      push( @$summary, $report_name );

      if ( getValue( $def, "perform_recluster" ) ) {
        my $recluster_name = $cluster_task_name . "_recluster";
        $config->{$recluster_name} = {
          class                => "CQS::UniqueR",
          perform              => 1,
          target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $recluster_name,
          rtemplate            => "../scRNA/scRNA_func.r;../scRNA/seurat_recluster.r",
          parameterFile1_ref   => [ $seurat_name, ".final.rds" ],
          parameterFile2_ref => [ $cluster_task_name, $cluster_file ],
          parameterSampleFile1 => $def->{recluster},
          parameterSampleFile2 => {
            recluster_celltypes => getValue( $def, "recluster_celltypes", "" ),
            Mtpattern           => getValue( $def, "Mtpattern" ),
            rRNApattern         => getValue( $def, "rRNApattern" ),
            Remove_Mt_rRNA      => getValue( $def, "Remove_Mt_rRNA" ),
            resolution          => getValue( $def, "resolution" ),
            pca_dims            => getValue( $def, "pca_dims" ),
            species             => getValue( $def, "species" ),
            markers_file        => getValue( $def, "markers_file" ),
            by_integration      => getValue( $def, "by_integration" ),
            by_sctransform      => getValue( $def, "by_sctransform" ),
            prefix              => $project_name,
            celltype_name     => $celltype_name,
            cluster_name      => $cluster_name
          },
          output_file_ext => ".recluster.rds",
          sh_direct       => 1,
          pbs             => {
            "email"     => $def->{email},
            "emailType" => $def->{emailType},
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "40gb"
          },
        };
        push( @$summary, $recluster_name );
      }

      if ( $def->{perform_antibody_vis} ) {
        addAntibodyTask( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name );
      }

      if(defined $clonotype_4_convert){
        my $clonotype_vis = $clonotype_4_convert;
        $clonotype_vis =~ s/3_convert/4_vis/ig;
        $clonotype_vis =~ s/4_convert/5_vis/ig;
        $config->{$clonotype_vis} = {
          class                      => "CQS::UniqueR",
          perform                    => 1,
          target_dir                 => $target_dir . "/" . $clonotype_vis,
          rtemplate                  => "../scRNA/scRNA_func.r;../scRNA/clonotype_vis.r",
          output_to_result_directory => 1,
          output_file_ext            => ".clonotype_vis.csv",
          parameterFile1_ref         => [ $seurat_name, ".final.rds" ],
          parameterFile2_ref         => [ $cluster_task_name, $cluster_file ],
          parameterFile3_ref         => [ $clonotype_4_convert ],
          sh_direct                  => 1,
          pbs                        => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "23",
            "mem"       => "10gb"
          },
        }
      }

      if ( $def->{perform_marker_dotplot} ) {
        my $biomarker_dotplot_task  = $cluster_task_name . "_biomarker_dotplot";
        my $marker_dotplot_clusters = getValue( $def, "marker_dotplot_clusters", {"all" => ["all"]});
        $config->{$biomarker_dotplot_task} = {
          class              => "CQS::UniqueR",
          perform            => 1,
          target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $biomarker_dotplot_task,
          rtemplate          => "../scRNA/scRNA_func.r;../scRNA/biomarker_dotplot.r",
          source             => $marker_dotplot_clusters,
          parameterFile1_ref => [ $seurat_name, ".final.rds" ],
          parameterFile2_ref => [ $cluster_task_name, $cluster_file ],
          parameterSampleFile1 => $marker_dotplot_clusters,
          parameterSampleFile2 => {
            cluster_name => getValue( $def, "marker_dotplot_clusters_name", "seurat_clusters" ),
            display_cluster_name => getValue( $def, "marker_dotplot_clusters_display_name", $cluster_name ),
            gene_number => getValue( $def, "marker_dotplot_gene_number", 20 ),
          },
          output_file => "parameterSampleFile1",
          output_file_ext    => ".count.files.csv",
          sh_direct          => 1,
          pbs                => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        };
        push( @$summary, $biomarker_dotplot_task );
      }

      if($def->{perform_curated_gene_dotplot}){
        my $curated_gene_dotplot_task  = $cluster_task_name . "_curated_gene_dotplot";
        my $curated_gene_def = $def->{curated_gene_dotplot};
        my ($expanded_gene_def, $clusters, $genes) = parse_curated_genes($curated_gene_def);

        $config->{$curated_gene_dotplot_task} = {
          class              => "CQS::UniqueR",
          perform            => 1,
          target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $curated_gene_dotplot_task,
          rtemplate          => "../scRNA/scRNA_func.r;../scRNA/curated_gene_dotplot.r",
          source             => $expanded_gene_def,
          parameterFile1_ref => [ $seurat_name, ".final.rds" ],
          parameterFile2_ref => [ $cluster_task_name, $cluster_file ],
          parameterSampleFile1 => $genes,
          parameterSampleFile2 => $clusters,
          parameterSampleFile3 => {
            cluster_name => "seurat_clusters",
            display_cluster_name => $cluster_name,
            by_sctransform => getValue($def, "by_sctransform"),
          },
          output_file => "parameterSampleFile1",
          output_file_ext    => ".count.files.csv",
          sh_direct          => 1,
          pbs                => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        };
        push( @$summary, $curated_gene_dotplot_task );
      }

      if ( $def->{t_cell_clusters} ) {
        my $tcell_clusters = $def->{t_cell_clusters};
        my $tcellTaskname  = $cluster_task_name . "_tcells";
        $config->{$tcellTaskname} = {
          class              => "CQS::UniqueR",
          perform            => 1,
          target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $tcellTaskname,
          rtemplate          => "../scRNA/extractTcells.r",
          parameterFile1_ref => [ $seurat_name, ".final.rds" ],
          parameterFile2     => $def->{marker_genes_file},
          parameterFile3_ref => [ $cluster_task_name, $cluster_file ],
          output_file_ext    => ".count.files.csv",
          rCode              => "celltype_name='" . $celltype_name . "'; cluster_name='" . join( ',', @$tcell_clusters ) . "'",
          sh_direct          => 1,
          pbs                => {
            "email"     => $def->{email},
            "emailType" => $def->{emailType},
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        };
        push( @$summary, $tcellTaskname );
      }

      my $perform_comparison = getValue( $def, "perform_comparison", 0 ) | getValue( $def, "perform_edgeR" );
      my $DE_by_sample = getValue( $def, "DE_by_sample" );
      my $DE_by_cell = getValue( $def, "DE_by_cell" );

      my @deByOptions = ();
      if ( getValue( $def, "DE_by_celltype" ) ) {
        push( @deByOptions, "DE_by_celltype" );
      }
      if ( getValue( $def, "DE_by_cluster" ) ) {
        push( @deByOptions, "DE_by_cluster" );
      }

      print("Perform_comparison=" . $perform_comparison , "\n");
      print("DE_by_sample=" . $DE_by_sample , "\n");
      
      if ( $perform_comparison & $DE_by_sample ) {
        for my $deByOption (@deByOptions) {
          my $DE_by_celltype = $deByOption eq "DE_by_celltype";
          addDeseq2BySampleTask( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name, 0, $DE_by_celltype, 0 );
        }
      }

      if ( $perform_comparison & $DE_by_cell ) {
        if ( defined $def->{"DE_cluster_pairs"} ) {
          addEdgeRTask( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name, 1, 0, 0 );
        }

        for my $deByOption (@deByOptions) {
          my $DE_by_celltype = $deByOption eq "DE_by_celltype";
          addEdgeRTask( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name, 0, $DE_by_celltype, 1 );
        }
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
