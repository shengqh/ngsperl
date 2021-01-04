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
  initDefaultValue( $def, "nFeature_cutoff_max", 5000 );
  initDefaultValue( $def, "nCount_cutoff",       500 );
  initDefaultValue( $def, "mt_cutoff",           20 );
  initDefaultValue( $def, "resolution",          0.8 );
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

sub addEdgeRTask {
  my ( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name, $bBetweenCluster, $DE_by_celltype, $DE_by_cell ) = @_;
  my $rCode = "pvalue=" . getValue( $def, "DE_pvalue" ) . ";useRawPvalue=" . getValue( $def, "DE_use_raw_pvalue" ) . ";foldChange=" . getValue( $def, "DE_fold_change" );

  my $edgeRtaskname  = $cluster_task_name . "_edgeR";
  my $groups         = undef;
  my $pairs          = undef;
  my $curClusterName = undef;
  if ($bBetweenCluster) {
    $edgeRtaskname  = $edgeRtaskname . "_betweenCluster_byCell";
    $curClusterName = getValue( $def, "DE_cluster_name" );
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
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
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
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push( @$summary, $vistaskname );

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
        "email"     => $def->{email},
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
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "23",
        "mem"       => "10gb"
      },
    };
    push( @$summary, $gseaTaskName );
  }

  return ($edgeRtaskname);
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

  if (defined $def->{vdj_files}){
    $config->{vdj_files} = $def->{vdj_files};
    addClonotypeMerge($config, $def, $summary, $target_dir, "clonotype_merge", ["vdj_files", "all_contig_annotations.json"]);
    addEnclone($config, $def, $summary, "clonotype_merge_enclone", $target_dir, ["clonotype_merge", ".json\$"] );
    addEncloneToClonotype($config, $def, $summary, $target_dir, "clonotype_enclone_to_clonotypes", "clonotype_merge_enclone", ["clonotype_merge", ".cdr3\$"]);
  }

  if (defined $def->{files}){
    my $split_hto_ref = undef;
    my $hto_sample_file = undef;
    if((defined $def->{perform_split_hto_samples}) and $def->{perform_split_hto_samples}) {
      my $r_script = undef;
      my $folder = undef;
      if ((defined $def->{split_hto_samples_by_cutoff}) and $def->{split_hto_samples_by_cutoff}) {
        $r_script = "../scRNA/split_samples_cutoff.r";
        $folder = "split_hto_samples_cutoff";
      } else {
        $r_script = "../scRNA/split_samples.r";
        $folder = "split_hto_samples_HTODemux";
      }
      $config->{"split_hto_samples"} = {
        class => "CQS::ProgramWrapperOneToOne",
        target_dir => "${target_dir}/$folder",
        interpretor => "R --vanilla -f ",
        program => $r_script,
        check_program => 1,
        option => "--args __FILE__ __OUTPUT__ " . getValue($def, "hto_regex", ""),
        source_arg => "",
        source_ref => [ "files" ],
        output_arg => "",
        output_file_prefix => ".HTO",
        output_file_ext => ".HTO.csv",
        output_to_same_folder => 0,
        can_result_be_empty_file => 0,
        sh_direct   => 1,
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$individual, "split_hto_samples" );

      $split_hto_ref = ["split_hto_samples", ".HTO.csv" ];

      $config->{"split_hto_samples_summary"} = {
        class => "CQS::ProgramWrapper",
        target_dir => "${target_dir}/${folder}_summary",
        interpretor => "R --vanilla -f ",
        program => "../scRNA/split_samples_summary.r",
        check_program => 1,
        option => "--args __FILE__ __OUTPUT__",
        source_arg => "",
        source_ref => $split_hto_ref,
        output_arg => "",
        output_file_prefix => ".HTO.summary",
        output_file_ext => ".HTO.summary.csv",
        sh_direct   => 1,
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      },

      push( @$individual, "split_hto_samples_summary" );

      if(defined $def->{HTO_samples}){
        $hto_sample_file = write_HTO_sample_file($def);
      }

      if(defined $def->{bam_files}){
        if (not defined $def->{HTO_samples}) {
          die "Define HTO_samples for split bam files";
        }

        $config->{HTO_samples} = $def->{HTO_samples};
        $config->{bam_files} = $def->{bam_files};
        $config->{"split_hto_bam"} = {
          class => "CQS::ProgramWrapperOneToOne",
          target_dir => "${target_dir}/split_hto_bam",
          interpretor => "python3",
          program => "../scRNA/split_samples.py",
          check_program => 1,
          option => "-o .",
          source_arg => "-i",
          source_ref => $split_hto_ref,
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
        push( @$individual, "split_hto_bam" );
      }

      if(defined $def->{vdj_json_files}){
        if (not defined $def->{HTO_samples}) {
          die "Define HTO_samples for split vdj json files";
        }

        $config->{HTO_samples} = $def->{HTO_samples};
        $config->{vdj_json_files} = $def->{vdj_json_files};
        $config->{"split_hto_clonotype"} = {
          class => "CQS::ProgramWrapperOneToManyFile",
          target_dir => "${target_dir}/split_hto_clonotype",
          interpretor => "python3",
          program => "../scRNA/clonotype_split.py",
          check_program => 1,
          option => "-o .",
          source_arg => "-i",
          source_ref => "vdj_json_files",
          parameterSampleFile2_arg => "-c",
          parameterSampleFile2_ref => $split_hto_ref,
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
        push( @$individual, "split_hto_clonotype" );

        addClonotypeMerge($config, $def, $summary, $target_dir, "split_hto_clonotype_merge", ["split_hto_clonotype", "all_contig_annotations.json"]);
        addEnclone($config, $def, $summary, "split_hto_clonotype_merge_enclone", $target_dir, ["split_hto_clonotype_merge", ".json\$"] );
        addEncloneToClonotype($config, $def, $summary, $target_dir, "split_hto_clonotype_merge_enclone_clonotypes", "split_hto_clonotype_merge_enclone", ["split_hto_clonotype_merge", ".cdr3\$"]);
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

    if ( getValue( $def, "perform_seurat" ) ) {
      my $additional_rmd_files = "Functions.Rmd";

      my $seurat_name = "seurat" . ( getValue( $def, "by_sctransform" ) ? "_sct" : "" ) . ( getValue( $def, "by_integration" ) ? "_igr" : "" ) . ( getValue( $def, "pool_sample" ) ? "_pool" : "" );
      $config->{$seurat_name} = {
        class                    => "CQS::UniqueRmd",
        perform                  => 1,
        target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $seurat_name,
        report_rmd_file          => "../scRNA/analysis.rmd",
        additional_rmd_files     => $additional_rmd_files,
        parameterSampleFile1_ref => "files",
        parameterSampleFile2     => {
          Mtpattern             => getValue( $def, "Mtpattern" ),
          rRNApattern           => getValue( $def, "rRNApattern" ),
          Remove_Mt_rRNA        => getValue( $def, "Remove_Mt_rRNA" ),
          nFeature_cutoff_min   => getValue( $def, "nFeature_cutoff_min" ),
          nFeature_cutoff_max   => getValue( $def, "nFeature_cutoff_max" ),
          nCount_cutoff         => getValue( $def, "nCount_cutoff" ),
          mt_cutoff             => getValue( $def, "mt_cutoff" ),
          resolution            => getValue( $def, "resolution" ),
          pca_dims              => getValue( $def, "pca_dims" ),
          species               => getValue( $def, "species" ),
          markers_file          => getValue( $def, "markers_file" ),
          details_rmd           => getValue( $def, "details_rmd" ),
          by_integration        => getValue( $def, "by_integration" ),
          by_sctransform        => getValue( $def, "by_sctransform" ),
          pool_sample           => getValue( $def, "pool_sample" ),
          batch_for_integration => getValue( $def, "batch_for_integration" ),
          hto_sample_file       => $hto_sample_file,
          prefix                => $taskName,
        },
        parameterSampleFile3 => $def->{"batch_for_integration_groups"},
        parameterSampleFile4 => $def->{"pool_sample_groups"},
        parameterSampleFile5_ref => $split_hto_ref,
        output_file_ext      => ".final.rds;.cluster.csv;.allmarkers.csv;.top10markers.csv;_ur.html",
        sh_direct            => 1,
        pbs                  => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary, $seurat_name );

      my $cluster_task_name = $seurat_name;
      my $cluster_file      = ".cluster.csv";
      my $celltype_name     = "cellactivity_clusters";
      my $cluster_name      = "seurat_cellactivity_clusters";

      addGeneTask( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name );

      if ( getValue( $def, "perform_rename_cluster" ) ) {
        my $rename_cluster = $seurat_name . "_rename_cluster";
        $config->{$rename_cluster} = {
          class                => "CQS::UniqueR",
          perform              => 1,
          target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $rename_cluster,
          rtemplate            => "../scRNA/renameCluster.r",
          parameterFile1_ref   => [ $seurat_name, ".final.rds" ],
          parameterSampleFile2 => $def->{rename_cluster},
          output_file_ext      => ".rename_cluster.csv;.rename_cluster.png",
          sh_direct            => 1,
          pbs                  => {
            "email"     => $def->{email},
            "emailType" => $def->{emailType},
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        };
        push( @$summary, $rename_cluster );

        $cluster_task_name = $rename_cluster;
        $cluster_file      = ".rename_cluster.csv";
        $celltype_name     = "renamed_cellactivity_clusters";
        $cluster_name      = "renamed_seurat_cellactivity_clusters";

        addGeneTask( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name );
      }

      if ( getValue( $def, "perform_recluster" ) ) {
        my $recluster_name = $seurat_name . "_recluster";
        $config->{$recluster_name} = {
          class                => "CQS::UniqueR",
          perform              => 1,
          target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $recluster_name,
          rtemplate            => "../scRNA/scRNArecluster.r",
          parameterFile1_ref   => [ $seurat_name, ".final.rds" ],
          parameterSampleFile2 => {
            recluster_celltypes => getValue( $def, "recluster_celltypes" ),
            Mtpattern           => getValue( $def, "Mtpattern" ),
            rRNApattern         => getValue( $def, "rRNApattern" ),
            Remove_Mt_rRNA      => getValue( $def, "Remove_Mt_rRNA" ),
            resolution          => getValue( $def, "resolution" ),
            pca_dims            => getValue( $def, "pca_dims" ),
            species             => getValue( $def, "species" ),
            markers_file        => getValue( $def, "markers_file" ),
            by_integration      => getValue( $def, "by_integration" ),
            by_sctransform      => getValue( $def, "by_sctransform" ),
            prefix              => $taskName,
          },
          output_file_ext => ".recluster.rds",
          sh_direct       => 1,
          pbs             => {
            "email"     => $def->{email},
            "emailType" => $def->{emailType},
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        };
        push( @$summary, $recluster_name );
      }

      if ( $def->{perform_antibody_vis} ) {
        addAntibodyTask( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name );
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

      my @deByOptions = ();
      if ( getValue( $def, "DE_by_celltype" ) ) {
        push( @deByOptions, "DE_by_celltype" );
      }
      if ( getValue( $def, "DE_by_cluster" ) ) {
        push( @deByOptions, "DE_by_cluster" );
      }

      my @deByOptions2 = ();
      if ( getValue( $def, "DE_by_cell" ) ) {
        push( @deByOptions2, "DE_by_cell" );
      }
      if ( getValue( $def, "DE_by_sample" ) ) {
        push( @deByOptions2, "DE_by_sample" );
      }

      if ( getValue( $def, "perform_edgeR" ) ) {
        if ( defined $def->{"DE_cluster_pairs"} ) {
          addEdgeRTask( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name, 1, 0, 0 );
        }

        for my $deByOption (@deByOptions) {
          my $DE_by_celltype = $deByOption eq "DE_by_celltype";

          for my $deByOption2 (@deByOptions2) {
            my $DE_by_cell = $deByOption2 eq "DE_by_cell";

            addEdgeRTask( $config, $def, $summary, $target_dir, $seurat_name, $cluster_task_name, $cluster_file, $celltype_name, $cluster_name, 0, $DE_by_celltype, $DE_by_cell );
          }
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
