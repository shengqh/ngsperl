#!/usr/bin/perl
package Pipeline::scRNASeq;

use strict;
use warnings;
use List::Util qw(first);
use File::Basename;
use Storable qw(dclone);
use File::Slurp;
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

our %EXPORT_TAGS = ( 'all' => [qw(initializeScRNASeqDefaultOptions 
  performScRNASeq 
  performScRNASeqTask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeScRNASeqDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );

  initDefaultValue( $def, "perform_scRNABatchQC", 1 );

  initDefaultValue( $def, "perform_individual_qc", 1 );

  initDefaultValue( $def, "perform_cellbender", 0 );

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
  initDefaultValue( $def, "DE_showVolcanoLegend", 0 );

  initDefaultValue( $def, "perform_clonotype_analysis", 0 );

  initDefaultValue( $def, "perform_SingleR", 0);
  initDefaultValue( $def, "perform_SignacX", 0);
  initDefaultValue( $def, "perform_Azimuth", 0);
  initDefaultValue( $def, "perform_DCATS", 0);

  initDefaultValue( $def, "perform_enclone_only",      0 );

  initDefaultValue( $def, "perform_seurat",      1 );
  
  initDefaultValue( $def, "seurat_walltime",          "24" );
  initDefaultValue( $def, "seurat_mem",          "120gb" );

  initDefaultValue( $def, "Mtpattern",           "^MT-|^Mt-" );
  initDefaultValue( $def, "rRNApattern",         "^Rp[sl][[:digit:]]|^RP[SL][[:digit:]]" );
  initDefaultValue( $def, "hemoglobinPattern",   "^HB[^(P)]|^Hb[^(p)]" );
  
  initDefaultValue( $def, "Remove_rRNA",         0 );
  initDefaultValue( $def, "Remove_MtRNA",        0 );
  initDefaultValue( $def, "regress_by_percent_mt", 1 );
  initDefaultValue( $def, "Remove_hemoglobin",   0 );
  
  initDefaultValue( $def, "nFeature_cutoff_min", 300 );
  initDefaultValue( $def, "nFeature_cutoff_max", 10000 );
  initDefaultValue( $def, "nCount_cutoff",       500 );
  initDefaultValue( $def, "mt_cutoff",           20 );
  initDefaultValue( $def, "resolution",          0.5 );
  initDefaultValue( $def, "details_rmd",         "" );

  initDefaultValue( $def, "by_sctransform", 1 );
  initDefaultValue( $def, "use_sctransform_v2", 1 );

  initDefaultValue( $def, "by_integration", 0 );
  if($def->{"by_integration"}){
    if(!defined $def->{integration_by_method}){
      if(getValue($def, "integration_by_fastmnn", 0)){
        $def->{integration_by_method} = "fastmnn";
      }elsif(getValue( $def, "integration_by_harmony", 1)){
        $def->{integration_by_method} = "harmony";
      }else{
        $def->{integration_by_method} = "seurat";
      }
    }

    $def->{integration_by_fastmnn} = $def->{integration_by_method} eq "fastmnn";
    $def->{integration_by_harmony} = $def->{integration_by_method} eq "harmony";
  }else{
    $def->{integration_by_fastmnn} = 0;
    $def->{integration_by_harmony} = 0;
  }

  my $pca_dims = $def->{by_sctransform}?30:20;
  #my $pca_dims = 50;
  initDefaultValue( $def, "pca_dims",            $pca_dims );

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

  if(defined $def->{"DE_by_sample"}){
    $def->{"DE_by_cell"} = $def->{"DE_by_sample"} ? 0: 1;
  }else{
    $def->{"DE_by_sample"} = 0;
    $def->{"DE_by_cell"} = 1;
  }

  initDefaultValue( $def, "DE_by_cell_filter_minTPM",         1 );
  initDefaultValue( $def, "DE_by_cell_filter_cellPercentage", 0.25 );

  initDefaultValue( $def, "DE_by_sample_min_cell_per_sample", 10 );
  
  initDefaultValue( $def, "perform_webgestalt", 0 );

  initDefaultValue( $def, "perform_CHETAH", 0 );

  initDefaultValue( $def, "perform_multires", 0 );
  
  initDefaultValue( $def, "perform_dynamic_cluster", 1 );
  initDefaultValue( $def, "dynamic_by_one_resolution", 0.2 );

  initDefaultValue( $def, "perform_sub_dynamic_cluster", 0);
  initDefaultValue( $def, "sub_dynamic_redo_harmony", 1);
  initDefaultValue( $def, "sub_dynamic_init_layer", "layer4");
  initDefaultValue( $def, "sub_dynamic_final_layer", "layer5");

  if(getValue($def, "species") ne "Hs"){
    if(getValue($def, "perform_SignacX", 0)){
      die "perform_SignacX should be 0 since the dataset is not from human species";
    }
  }

  initDefaultValue( $def, "perform_dynamic_cluster_signacX", getValue($def, "perform_SignacX", 0) );

  # initDefaultValue( $def, "dynamic_layer_umap_min_dist", {
  #   "layer0" => 0.3,
  #   "layer1" => 0.3,
  #   "layer2" => 0.2,
  #   "layer3" => 0.1,
  #   "layer4" => 0.05,
  # } );
  
  initDefaultValue( $def, "perform_fix_resolution", 0 );
  #initDefaultValue( $def, "remove_subtype", "T cells,Fibroblasts,Neurons,Macrophages,Dendritic cells"),
  initDefaultValue( $def, "remove_subtype", "T cells,B cells,Plasma cells,Fibroblasts,Neurons,Epithelial cells,Endothelial cells,Macrophages,Dendritic cells,Ciliated cells"),
  initDefaultValue( $def, "best_resolution_min_markers", 20);

  initDefaultValue( $def, "perform_decontX", 0);
  initDefaultValue( $def, "remove_decontX", 0);
  initDefaultValue( $def, "remove_decontX_by_contamination", 0);

  return $def;
}


sub getScRNASeqConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeScRNASeqDefaultOptions($def);

  my $project_name = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $tasks = [@$individual, @$summary];

  my $clono_key = "clono_index";

  my $target_dir      = $def->{target_dir};
  my $groups_ref      = defined $def->{groups} ? "groups" : undef;
  my $aligner         = $def->{aligner};
  my $star_option     = $def->{star_option};
  my $count_table_ref = "files";

  $config->{bam_files} = $def->{bam_files};

  my $essential_gene_task;
  
  if(getValue($def, "perform_essential_gene", 1)){
    $essential_gene_task = add_essential_gene($config, $def, $tasks, $target_dir);
  } 

  my $perform_comparison = getValue( $def, "perform_comparison", 0 );
  if(getValue( $def, "perform_edgeR" )){
    $perform_comparison = 1;
  }
  my $DE_by_sample = getValue( $def, "DE_by_sample" );
  my $DE_by_cell = getValue( $def, "DE_by_cell" );

  my @deByOptions = ();
  if ( getValue( $def, "DE_by_celltype" ) ) {
    push( @deByOptions, "DE_by_celltype" );
  }
  if ( getValue( $def, "DE_by_cluster" ) ) {
    push( @deByOptions, "DE_by_cluster" );
  }

  print("perform_comparison=" . $perform_comparison , "\n");
  print("DE_by_sample=" . $DE_by_sample , "\n");

  my $perform_split_hto_samples = getValue($def, "perform_split_hto_samples", 0);
  if($perform_split_hto_samples){
    my $hto_samples = getValue($def, "HTO_samples");
    $def->{hto_file_names} = [sort keys %$hto_samples];
  }
  my $hto_file_names = $def->{hto_file_names};

  my $perform_clonotype_analysis = getValue($def, "perform_clonotype_analysis", 0);
  my $clonotype_ref = undef;

  my $perform_arcasHLA = getValue($def, "perform_arcasHLA", 0);
  my $perform_strelka2 = getValue($def, "perform_strelka2", 0);

  my $bam_ref = undef;
  my $hla_merge = undef;
  my $filter_config_file = getValue($def, "filter_config_file", "");
  if($filter_config_file eq ""){
    $filter_config_file = undef;
  }
  my $sctk_ref = undef;
  my $signacX_ref = undef;
  my $singleR_ref = undef;
  my $azimuth_ref = undef;
  my $decontX_ref = undef;
  my $doublet_finder_ref = undef;

  if (defined $def->{files}){
    if(defined $def->{atac_files}){
      $config->{atac_files} = $def->{atac_files};
      my $fragment_cells_task = undef;
      if(getValue($def, "perform_fragment_cells", 1)){
        $fragment_cells_task = "fragment_cells";
        add_fragment_cells($config, $def, $tasks, $target_dir, $fragment_cells_task, "atac_files");
      }

      my $multiome_task = "multiome_qc";
      add_multiome_qc($config, $def, $tasks, $target_dir, $multiome_task, undef, $fragment_cells_task);
    }

    my @report_files = ();
    my @report_names = ();
    my $hto_task = undef;
    my $hto_ref = undef;
    my $hto_sample_file = undef;
    my $hto_summary_task = undef;
    my $raw_files_def = "raw_files";
    my $files_def = "files";
    my $filtered_files_def = "files";

    my $perform_cellbender = getValue($def, "perform_cellbender", 0);

    my $perform_decontX = getValue($def, "perform_decontX", 0);
    my $remove_decontX = $perform_decontX && getValue($def, "remove_decontX", 0);
    my $remove_decontX_by_contamination = getValue($def, "remove_decontX_by_contamination", 0.25);

    my $prefix = "";

    my $raw_individual_qc_task = undef;
    my $qc_report_task = undef;

    my $perform_sctk = getValue($def, "perform_sctk", 0);
    my $remove_doublets = getValue($def, "remove_doublets", 0);
    my $perform_individual_qc = getValue($def, "perform_individual_qc", 1);

    my $decontX_task = undef;
    my $decontX_counts_ref = undef;
    if ($perform_decontX){
      $decontX_task = "decontX";
      add_decontX($config, $def, $tasks, $target_dir, $decontX_task, $filtered_files_def, $raw_files_def, {}, 1);
      $decontX_ref = [$decontX_task, ".meta.rds"];
      $decontX_counts_ref = [$decontX_task, ".counts.rds"];
      if($remove_decontX){
        $files_def = $decontX_counts_ref;
        $prefix = "decontX_";
      }
    }

    if($perform_cellbender){
      my $cellbender_prefix = "cellbender";
      $files_def = add_cellbender_v2($config, $def, $tasks, $target_dir, $cellbender_prefix, $filtered_files_def, $raw_files_def, $decontX_counts_ref );

      if($remove_decontX){
        $prefix = $prefix . "cellbender_";
      }else{
        $prefix = "cellbender_";
      }

      $remove_decontX = 0;
    }

    if ( $perform_sctk ){
      my $sctk_task = $prefix . "sctk";
      add_sctk($config, $def, $tasks, $target_dir, $sctk_task, $files_def);
      $sctk_ref = [$sctk_task, ".meta.rds"];
    }


    if($remove_doublets){
      my $doublets_ref = undef;
      my $doublet_column = undef;
      my $nodoublets_task = "${prefix}nd";
      if (defined $doublet_finder_ref){
        $doublets_ref = $doublet_finder_ref;
        $doublet_column = getValue($def, "remove_doublet_column", getValue($def, "doublet_column", "DF.classifications_highest"));
      }elsif ($perform_sctk){
        $doublets_ref = $sctk_ref;
        $doublet_column = getValue($def, "remove_doublet_column", getValue($def, "doublet_column", "doubletFinder_doublet_label_resolution_1.5"));
      }

      if(defined $doublets_ref){  
        add_remove_doublets($config, $def, $tasks, $target_dir, $nodoublets_task, $files_def, $doublets_ref, $doublet_column);
        
        $files_def = [$nodoublets_task, ".counts.rds"];
        $prefix = "${prefix}nd_";
      }
    }

    if ( $perform_individual_qc ){
      ($raw_individual_qc_task, $qc_report_task, $signacX_ref, $singleR_ref, $azimuth_ref) = add_individual_qc_tasks($config, $def, $tasks, $target_dir, $project_name, $prefix, $filter_config_file, $files_def, $decontX_ref, $sctk_ref);
    }

    if( $def->{"perform_individual_dynamic_qc"} ){
      my $sct_str = get_sct_str($def);
      my $raw_individual_dynamic_qc_task = "${prefix}raw_dynamic_qc${sct_str}";
      add_individual_dynamic_qc($config, $def, $tasks, $target_dir, $raw_individual_dynamic_qc_task, $filter_config_file, $files_def, $essential_gene_task);
    }

    my $files = $def->{files};
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

    my $hto_raw_file_ref = undef;
    my $raw_files = $def->{raw_files};
    if(defined $raw_files){
      $hto_raw_file_ref = "raw_files";
      if(defined $hto_file_names){
        my $hto_raw_files = {};
        for my $hto_name (@$hto_file_names){
          $hto_raw_files->{$hto_name} = $raw_files->{$hto_name};
        }
        $config->{hto_raw_files} = $hto_raw_files;
        $hto_raw_file_ref = "hto_raw_files";
      }
    }

    # if(getValue($def, "perform_scDblFinder", 0)){
    #   add_scDblFinder($config, $def, $tasks, $target_dir, "scDblFinder", "files" );
    #   #$files_def = "scDblFinder";
    # }

    if( $perform_split_hto_samples ) {
      my $preparation_task = add_hto_samples_preparation($config, $def, $tasks, $target_dir, $hto_file_ref, $hto_raw_file_ref);
      $hto_file_ref = [ $preparation_task, ".hto.rds"];

      if(defined $def->{HTO_samples}){
        $hto_sample_file = write_HTO_sample_file($def);
      }

      if(getValue($def, "split_hto_samples_by_GMM_demux", 0)){
        my $gmm_demux_task = add_hto_gmm_demux($config, $def, $tasks, $target_dir, $hto_file_ref, $hto_sample_file);
        #my $hto_summary_task = add_hto_summary($config, $def, $tasks, $target_dir, [$gmm_demux_task, ".HTO.csv"]);
      }

      my $hto_task = add_hto($config, $def, $tasks, $target_dir, $hto_file_ref, $hto_sample_file);
      $hto_ref = [ $hto_task, ".HTO.csv" ];

      my $hto_bam_ref = $hto_ref;

      #my $hto_summary_task = add_hto_summary($config, $def, $tasks, $target_dir, $hto_ref);

      #push (@report_files, ($hto_summary_task, ".HTO.summary.global.png"));
      #push (@report_names, "hto_summary_png");

      my $perform_souporcell = getValue($def, "perform_souporcell", 0);
      if($perform_souporcell) {
        my $hto_souporcell_task = add_souporcell($config, $def, $tasks, $target_dir, $preparation_task);
        my $hto_integration_task = add_souporcell_integration($config, $def, $tasks, $target_dir, $hto_souporcell_task, $hto_ref);
        $hto_ref = [ $hto_integration_task, ".meta.rds" ];
        $hto_bam_ref = [ $hto_integration_task, ".HTO.csv" ];
      }

      if($perform_arcasHLA || $perform_strelka2){
        $bam_ref = add_hto_bam($config, $def, $tasks, $target_dir, $hto_bam_ref);
      }

      if($perform_arcasHLA){
        $hla_merge = addArcasHLA($config, $def, $tasks, $target_dir, $project_name, "", $bam_ref);        
      }

      if($perform_clonotype_analysis){
        my $split_task = add_clonotype_split($config, $def, $tasks, $target_dir, $hto_bam_ref, $clono_key);
        $clonotype_ref = [$split_task, "all_contig_annotations.json"];
      }
    }else{
      if($perform_arcasHLA || $perform_strelka2){
        getValue($def, "bam_files");
        $bam_ref = "bam_files";
      }
    
      if($perform_arcasHLA){
        if (defined $def->{singleend_bam_files}){
          $config->{singleend_bam_files} = $def->{singleend_bam_files};
          $hla_merge = addArcasHLA($config, $def, $tasks, $target_dir, $project_name, "", $bam_ref, "singleend_bam_files");        
        }else{
          $hla_merge = addArcasHLA($config, $def, $tasks, $target_dir, $project_name, "", $bam_ref);        
        }
      }
    }

    if($perform_strelka2){
      my $strelka2_task = "strelka2";
      add_strelka2($config, $def, $tasks, $target_dir, $strelka2_task, $bam_ref);
    }

    my $merge_task = undef;
    my $enclone_task = undef;
    my $clonotype_convert = undef;
    my $clonotype_db = undef;
    if ($perform_clonotype_analysis){
      if ((not defined $def->{files}) || (not $perform_split_hto_samples)) {
        $config->{vdj_json_files} = getValue($def, "vdj_json_files");
        $clonotype_ref = ["vdj_json_files", "all_contig_annotations.json"];
      }
      
      $merge_task = addClonotypeMerge($config, $def, $tasks, $target_dir, "clonotype". get_next_index($def, $clono_key) . "_merge", $clonotype_ref);
      $enclone_task = addEnclone($config, $def, $tasks, "clonotype". get_next_index($def, $clono_key) . "_enclone", $target_dir, [$merge_task, ".json\$"] );
      $clonotype_convert = addEncloneToClonotype($config, $def, $tasks, $target_dir, "clonotype". get_next_index($def, $clono_key) . "_convert", [$enclone_task, ".pchain4.csv"], [$merge_task, ".cdr3\$"]);
      $clonotype_db = addClonotypeDB($config, $def, $tasks, $target_dir, "clonotype". get_next_index($def, $clono_key) . "_db", $clonotype_convert);
      my $clonotype_consensus = addEncloneToConsensus($config, $def, $tasks, $target_dir, "clonotype". get_next_index($def, $clono_key) . "_consensus", [$enclone_task, ".pchain4.pcell.csv"], [$merge_task, ".cdr3\$"]);
      my $immunarch_task = addConsensusToImmunarch($config, $def, $tasks, $target_dir, "clonotype". get_next_index($def, $clono_key) . "_immunarch", $clonotype_consensus);
    }

    if ( getValue( $def, "perform_scRNABatchQC" ) ) {
      add_scRNABatchQC($config, $def, $tasks, $target_dir);
    }

    my $seurat_rawdata;
    my $is_preprocessed;
    if ( getValue( $def, "perform_seurat" ) ) {
      my $seurat_task = undef;
      my $reduction = undef;
      my $obj_ref = undef;
      my $localization_ref = undef;

      if(getValue($def, "rawdata_from_object", 0)){
        $seurat_task = "files";
        $reduction = "pca";
        $obj_ref = $seurat_task;
        $localization_ref = $seurat_task;
      }else{
        if(getValue($def, "rawdata_from_qc", 0)){
          if(!defined $raw_individual_qc_task){
            die("trying to build rawdata from qc, please set perform_individual_qc => 1 in your configuration file");
          }
          my $sct_str = get_sct_str($def);
          $seurat_rawdata = "${prefix}seurat_rawdata_postqc${sct_str}";
          add_seurat_rawdata($config, $def, $tasks, $target_dir, $seurat_rawdata, $hto_ref, $hto_sample_file, $files_def, undef, undef, $raw_individual_qc_task, $decontX_ref );
          $is_preprocessed = 1;
        }else{
          $seurat_rawdata = "${prefix}seurat_rawdata";

          if (getValue($def, "merge_seurat_object", 0)){
            add_seurat_merge_object($config, $def, $tasks, $target_dir, $seurat_rawdata, $files_def, undef, undef, 0, {});
          }else{
            add_seurat_rawdata($config, $def, $tasks, $target_dir, $seurat_rawdata, $hto_ref, $hto_sample_file, $files_def, undef, undef );
          }
          $is_preprocessed = 0;
        }

        push (@report_files, ($seurat_rawdata, "rawobj.rds"));
        push (@report_names, "raw_obj");

        ($seurat_task, $reduction) = add_seurat($config, $def, $tasks, $target_dir, $seurat_rawdata, $essential_gene_task, 0, $is_preprocessed, $prefix, $filter_config_file);
        $obj_ref = [$seurat_task, ".final.rds"];

        push (@report_files, ($seurat_task, ".final.png", 
          $seurat_task, ".qc.1.png", 
          $seurat_task, ".qc.2.png", 
          $seurat_task, ".qc.3.png", 
          $seurat_task, ".qc.4.png", 
          $seurat_task, ".sample_cell.csv"));
        push (@report_names, ("seurat_merge_png", "seurat_qc_1_png", "seurat_qc_2_png", "seurat_qc_3_png", "seurat_qc_4_png", 
          "sample_cell_csv"));

        $localization_ref = [ $seurat_task, ".final.rds" ];
      }

      if(!defined $signacX_ref) {
        if (getValue( $def, "perform_SignacX", 0 ) ) {
          my $signacX_task = $seurat_task . "_SignacX";
          add_signacx( $config, $def, $tasks, $target_dir, $project_name, $signacX_task, $obj_ref, $reduction );
          $signacX_ref = [$signacX_task, ".meta.rds"];
        }
      }

      if(!defined $singleR_ref) {
        if (getValue( $def, "perform_SingleR", 0 ) ) {
          my $singleR_task = $seurat_task . "_SingleR";
          my $cur_options = {
            task_name => $def->{task_name},
            reduction => $reduction, 
          };
          add_singleR_cell( $config, $def, $tasks, $target_dir, $singleR_task, $obj_ref, $cur_options );
          $singleR_ref = [$singleR_task, ".meta.rds"];
        }
      }

      if(!defined $azimuth_ref){
        if (getValue( $def, "perform_Azimuth", 0 ) ) {
          my $azimuth_task = $seurat_task . "_Azimuth";
          my $cur_options = {
            task_name => $def->{task_name},
            reduction => $reduction, 
          };
          add_azimuth( $config, $def, $tasks, $target_dir, $azimuth_task, $obj_ref, $cur_options);
          $azimuth_ref = [ $azimuth_task, ".meta.rds" ];
        }
      }

      my $celltype_task = undef;
      my $celltype_name = undef;

      if(getValue($def, "perform_dynamic_cluster")){
        my $raw_dynamicKey = $seurat_task . "_dr" . getValue($def, "dynamic_by_one_resolution");

        my $dynamicKey = $raw_dynamicKey . "_";

        $def->{$dynamicKey} = 0;

        my $scDynamic_task;

        my $by_individual_sample=0;
        my $by_column=undef;
        my $by_harmony;
        
        if(getValue($def, "dynamic_by_harmony", 0)){
          $by_harmony = 1;
          $scDynamic_task = $dynamicKey . get_next_index($def, $dynamicKey) . "_call_harmony";
        } else {
          $by_harmony = 0;
          $scDynamic_task = $dynamicKey . get_next_index($def, $dynamicKey) . "_call";
        }
        addDynamicCluster($config, $def, $tasks, $target_dir, $scDynamic_task, $seurat_task, $essential_gene_task, $reduction, $by_individual_sample, $by_column, $by_harmony);

        my $meta_ref = [$scDynamic_task, ".meta.rds"];
        my $call_files_ref = [$scDynamic_task, ".iter_png.csv"];

        if (defined $sctk_ref or defined $signacX_ref or defined $singleR_ref){
          my $validation_task = $scDynamic_task . "_validation";
          add_celltype_validation( $config, $def, $tasks, $target_dir, $validation_task, $seurat_task, $meta_ref, $call_files_ref, "layer4", ".dynamic_call_validation.html", 0, $signacX_ref, $singleR_ref, $sctk_ref, $decontX_ref, $azimuth_ref);
        }

        if(defined $def->{bubble_files}){
          add_bubble_files($config, $def, $tasks, $target_dir, $scDynamic_task . "_bubble_files", $seurat_task, $meta_ref, "layer4", "layer4_clusters", ".dynamic_layer4_bubbles.html" );
        }

        if(defined $def->{bubble_plots}){
          #add_bubble_plots($config, $def, $tasks, $target_dir, $scDynamic_task . "_bubblemap_iter1", $seurat_task, $meta_ref, "iter1", "iter1_clusters", ".dynamic_iter1_dot.html" );
          add_bubble_plots($config, $def, $tasks, $target_dir, $scDynamic_task . "_bubblemap", $seurat_task, $meta_ref, "layer4", "layer4_clusters", ".dynamic_layer4_dot.html" );
        }

        if(getValue($def, "perform_individual_dynamic_cluster", 0)){
          my $individual_scDynamic_task = $raw_dynamicKey . "_individual";
          addDynamicCluster($config, $def, $tasks, $target_dir, $individual_scDynamic_task, $seurat_task, $essential_gene_task, "pca", 1);

          my $clustree_task = $individual_scDynamic_task . "_clustree";
          add_clustree_rmd($config, $def, $tasks, $target_dir, $clustree_task, $individual_scDynamic_task, $scDynamic_task);
        }

        if(getValue($def, "perform_sub_dynamic_cluster")){
          my $sub_scDynamic_task = $scDynamic_task . "_sub";
          addSubDynamicCluster($config, $def, $tasks, $target_dir, $sub_scDynamic_task, $seurat_task, $meta_ref, $essential_gene_task, $reduction, 0);
        }


        # if (getValue( $def, "perform_SingleR", 0 ) ) {
        #   my $singleR_task = $scDynamic_task . "_SingleR";
        #   my $cur_options = {
        #     task_name => $def->{task_name},
        #     reduction => $reduction, 
        #     celltype_layer => "layer4",
        #     celltype_cluster => "layer4_clusters"
        #   };
        #   add_singleR( $config, $def, $tasks, $target_dir, $singleR_task, $obj_ref, $meta_ref, $cur_options );
        # }

        if(getValue($def, "perform_dynamic_subcluster")){
          my $subcluster_task = $dynamicKey . get_next_index($def, $dynamicKey) . "_subcluster" . (getValue($def, "subcluster_redo_harmony") ? "_rh" : "") ;
          
          my $cur_options = {
            reduction => $reduction, 
            celltype_layer => getValue($def, "dynamic_subcluster_init_celltype_layer", "layer4"),
            celltype_cluster => getValue($def, "dynamic_subcluster_init_celltype_cluster", "layer4_clusters")
          };

          my $rename_map = $def->{"dynamic_rename_map"};
          if(defined $def->{dynamic_delete_celltypes}){
            my $delete_cts = $def->{dynamic_delete_celltypes};
            for my $dct (@$delete_cts){
              if(!defined $rename_map){
                $rename_map = {};
              }
              $rename_map->{$dct} = {
                from => $dct,
                cluster => -1,
                to => "DELETE",
              };
            }
          }

          $subcluster_task = addSubCluster($config, $def, $tasks, $target_dir, $subcluster_task, $obj_ref, $meta_ref, $essential_gene_task, $cur_options, $rename_map, ".dynamic_subcluster.html", $signacX_ref, $singleR_ref,  $azimuth_ref);
          $meta_ref = [$subcluster_task, ".meta.rds"];

          if(getValue($def, "perform_dynamic_choose")) {
            my $choose_task = $dynamicKey . get_next_index($def, $dynamicKey) . "_choose";
            my $table = getValue($def, "dynamic_subclusters_table");
            addSubClusterChoose($config, $def, $tasks, $target_dir, $choose_task, $obj_ref, $meta_ref, $subcluster_task, $essential_gene_task, $cur_options, $table, ".dynamic_choose.html");
            $obj_ref = [ $choose_task, ".final.rds" ];
            $meta_ref = [ $choose_task, ".meta.rds" ];

            if (defined $sctk_ref or defined $signacX_ref or defined $singleR_ref){
              my $validation_task = $choose_task . "_validation";
              add_celltype_validation( $config, $def, $tasks, $target_dir, $validation_task, $obj_ref, $meta_ref, undef, "seurat_cell_type", ".dynamic_choose_validation.html", 1, $signacX_ref, $singleR_ref, $sctk_ref, $decontX_ref, $azimuth_ref );
            }

            $celltype_task = $choose_task;

            if($def->{perform_DCATS}){
              my $dcats_task = $choose_task . "_DCATS";
              add_dcats($config, $def, $summary, $target_dir, $dcats_task, [$celltype_task, ".final.rds"]);
            }

            if(defined $def->{bubble_files}){
              add_bubble_files($config, $def, $tasks, $target_dir, $choose_task . "_bubble_files", $choose_task, undef, undef, undef, ".dynamic_choose_bubbles.html" );
            }

            if(defined $def->{bubble_plots}){
              my $bubble_task = $choose_task . "_bubblemap";
              add_bubble_plots($config, $def, $tasks, $target_dir, $bubble_task, $choose_task, undef, undef, undef, ".dynamic_choose_dot.html");
            }

            if(defined $clonotype_convert) {
              if( getValue($def, "perform_gliph2", 0) ) {
                my $gliph2_task = add_gliph2($config, $def, $tasks, $target_dir, $meta_ref, $clonotype_convert, $hla_merge);
              }
              my $clonotype_consensus = addEncloneToConsensus($config, $def, $tasks, $target_dir, "clonotype". get_next_index($def, $clono_key) . "_consensus_tcell", [$enclone_task, ".pchain4.pcell.csv"], [$merge_task, ".cdr3\$"], [ $choose_task, ".meta.csv" ]);
              my $immunarch_task = addConsensusToImmunarch($config, $def, $tasks, $target_dir, "clonotype". get_next_index($def, $clono_key) . "_immunarch_tcell", $clonotype_consensus);
            }

            addComparison($config, $def, $tasks, $target_dir, $choose_task, $choose_task, "", "cell_type", "seurat_cell_type", "subumap");

            $localization_ref = $obj_ref;

            if(defined $def->{groups}){
              my $group_umap_task = $choose_task . "_group_umap";
              add_group_umap($config, $def, $tasks, $target_dir, $group_umap_task, [$choose_task, ".final.rds"]);
            }

            if(defined $clonotype_convert){
              addClonotypeCluster($config, $def, $tasks, $target_dir, $clonotype_convert . "_cluster", $clonotype_convert, $meta_ref, ".clonotype.cluster.csv,.clonotype.sub.cluster.csv");
              addClonotypeVis($config, $def, $tasks, $target_dir, $clonotype_convert . "_vis", $obj_ref, undef, $clonotype_convert);
            }

            if(defined $clonotype_db){
              addClonotypeCluster($config, $def, $tasks, $target_dir, $clonotype_db . "_cluster", $clonotype_db, $meta_ref, ".clonotype.db.cluster.csv,.clonotype.sub.db.cluster.csv");
            }

            # if(defined $df_task){
            #   my $doublet_check_task = $dynamicKey . get_next_index($def, $dynamicKey) . "_doublet_check";
            #   add_doublet_check($config, $def, $tasks, $target_dir, $doublet_check_task, $obj_ref, $df_task );
            # }
          }
        }
      }

      if(getValue($def, "perform_localization_genes_plot", 0)){
        my $gene_localization_map_task  = $localization_ref->[0] . "_gene_localization_map";

        my $all_genes=[];
        if(defined $def->{localization_genes} && $def->{localization_genes} ne ""){
          $all_genes = $def->{localization_genes};
        }

        if(defined $def->{localization_genes_file}){
          die "FILE_NOT_FOUND: localization_genes_file " . $def->{localization_genes_file} if (! -e $def->{localization_genes_file});
          my @genes_in_file = read_file($def->{localization_genes_file}, chomp => 1);
          for my $gene (@genes_in_file){
            push(@$all_genes, $gene);
          }
        }

        die "No localization_genes or localization_genes_file defined." if(scalar(@$all_genes) == 0);

        my $genes = {
          $def->{task_name} => $all_genes,
        };
        #print(Dumper($genes));
        $config->{$gene_localization_map_task} = {
          class              => "CQS::UniqueR",
          perform            => 1,
          target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $gene_localization_map_task,
          rtemplate          => "../scRNA/scRNA_func.r;../scRNA/gene_localization_map.r",
          source             => $genes,
          parameterFile1_ref => $localization_ref,
          parameterSampleFile1 => $genes,
          parameterSampleFile2 => $def->{groups},
          output_file_ext    => ".figure.files.csv",
          sh_direct          => 1,
          pbs                => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "12",
            "mem"       => "40gb"
          },
        };
        push( @$summary, $gene_localization_map_task );
      }

      if(getValue($def, "perform_localization_gene_ratio_plot", 0)){
        my $gene_ratio_localization_map_task  = $seurat_task . "_gene_ratio_localization_map";
        my $genes = {
          $def->{task_name} => getValue($def, "localization_gene_ratio"),
        };
        #print(Dumper($genes));
        $config->{$gene_ratio_localization_map_task} = {
          class              => "CQS::UniqueR",
          perform            => 1,
          target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $gene_ratio_localization_map_task,
          rtemplate          => "../scRNA/scRNA_func.r;../scRNA/gene_ratio_localization_map.r",
          source             => $genes,
          parameterFile1_ref => [ $seurat_task, ".final.rds" ],
          parameterSampleFile1 => $genes,
          parameterSampleFile2 => $def->{groups},
          output_file_ext    => ".figure.files.csv",
          sh_direct          => 1,
          pbs                => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "12",
            "mem"       => "40gb"
          },
        };
        push( @$summary, $gene_ratio_localization_map_task );
      }

      if(getValue($def, "perform_multires", 0)){
        my $multiresKey = $seurat_task . "_multires";
        $def->{$multiresKey} = 0;
        my $multires_task = $seurat_task . "_multires" . get_next_index($def, $multiresKey) . "_call";

        my $obj_ref = [$seurat_task, ".rds"];
        $config->{$multires_task} = {
          class                    => "CQS::UniqueR",
          perform                  => 1,
          target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $multires_task,
          rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_multires.r",
          rReportTemplate          => "../scRNA/seurat_multires.rmd;reportFunctions.R",
          parameterFile1_ref => $obj_ref,
          parameterSampleFile1    => {
            task_name             => getValue( $def, "task_name" ),
            Mtpattern             => getValue( $def, "Mtpattern" ),
            rRNApattern           => getValue( $def, "rRNApattern" ),
            Remove_rRNA           => getValue( $def, "Remove_rRNA" ),
            Remove_MtRNA          => getValue( $def, "Remove_MtRNA" ),
            pca_dims              => getValue( $def, "pca_dims" ),
            by_sctransform        => getValue( $def, "by_sctransform" ),
            by_integration        => getValue( $def, "by_integration" ),
            reduction             => $reduction,
            species               => getValue( $def, "species" ),
            db_markers_file       => getValue( $def, "markers_file" ),
            curated_markers_file  => getValue( $def, "curated_markers_file", "" ),
            annotate_tcell        => getValue( $def, "annotate_tcell", 0),
            remove_subtype        => getValue( $def, "remove_subtype", ""),
            HLA_panglao5_file     => getValue( $def, "HLA_panglao5_file", "" ),
            tcell_markers_file    => getValue( $def, "tcell_markers_file", ""),
            summary_layer_file    => $def->{summary_layer_file},
            bubblemap_file        => $def->{bubblemap_file},
            bubblemap_use_order   => getValue($def, "bubblemap_use_order", 0),
            plot_width            => getValue($def, "multires_plot_width", 9900),
            plot_height           => getValue($def, "multires_plot_height", 6000),
            layer                 => getValue($def, "multires_layer", "Layer4")
          },
          parameterSampleFile2    => getValue($def, "multires_combine_cell_types", {
            "NK/T cells" => ["NK cells", "T cells"],
            "Monocytes" => ["Macrophages", "Dendritic cells"],
          }),
          output_file_ext      => ".meta.rds",
          output_other_ext  => ".resolutions.csv",
          sh_direct            => 1,
          pbs                  => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "23",
            "mem"       => getValue($def, "seurat_mem")
          },
        };
        push( @$summary, $multires_task );
        my $meta_ref = [$multires_task, ".meta.rds"];

        if(getValue($def, "perform_multires_subcluster", 0)){
          my $subcluster_task = $seurat_task . "_multires" . get_next_index($def, $multiresKey) . "_subcluster";
          my $multires_resolution = getValue( $def, "multires_resolution" );

          my $celltype_cluster;
          if(getValue( $def, "by_sctransform" )){
            $celltype_cluster = "SCT_snn_res." . $multires_resolution;
          }else{
            $celltype_cluster = "RNA_snn_res." . $multires_resolution;
          }
          if(getValue( $def, "by_integration" ) & (! getValue($def, "integration_by_harmony"))) {
            $celltype_cluster = "integrated_snn_res." . $multires_resolution;
          }
          my $multires_celltype = $celltype_cluster . "_celltype_summary";

          if (defined $sctk_ref or defined $signacX_ref or defined $singleR_ref){
            my $validation_task = $multires_task . "_validation";
            add_celltype_validation( $config, $def, $tasks, $target_dir, $validation_task, $seurat_task, $meta_ref, undef, $celltype_cluster . "_celltype", ".multires_call_validation.html", 0, $signacX_ref, $singleR_ref, $sctk_ref, $decontX_ref, $azimuth_ref );
          }

          my $cur_options = {
            reduction => $reduction, 
            celltype_layer => $multires_celltype,
            celltype_cluster => $celltype_cluster,
          };
          my $rename_map = $def->{"multires_rename_map"};

          addSubCluster($config, $def, $tasks, $target_dir, $subcluster_task, $obj_ref, $meta_ref, $essential_gene_task, $cur_options, $rename_map, ".multires_subcluster.html", $signacX_ref, $singleR_ref, $azimuth_ref);
          $meta_ref = [$subcluster_task, ".meta.rds"];

          if(getValue($def, "perform_multires_choose", 0)) {
            my $choose_task = $seurat_task . "_multires" . get_next_index($def, $multiresKey) . "_choose";
            my $table = getValue($def, "multires_subclusters_table");
            addSubClusterChoose($config, $def, $tasks, $target_dir, $choose_task, $obj_ref, $meta_ref, $subcluster_task, $essential_gene_task, $cur_options, $table, ".multires_choose.html");

            if (defined $sctk_ref or defined $signacX_ref or defined $singleR_ref){
              my $validation_task = $choose_task . "_validation";
              add_celltype_validation( $config, $def, $tasks, $target_dir, $validation_task, [$choose_task, ".final.rds"], [$choose_task, "meta.rds"], undef, "seurat_cell_type", ".multires_choose_validation.html", 1, $signacX_ref, $singleR_ref, $sctk_ref, $decontX_ref, $azimuth_ref );
            }

            $celltype_task = $choose_task;

            my $obj_ref = [$choose_task, ".final.rds"];
            if(defined $def->{groups}){
              my $group_umap_task = $seurat_task . "_multires" . get_next_index($def, $multiresKey) . "_group_umap";
              add_group_umap($config, $def, $tasks, $target_dir, $group_umap_task, $obj_ref);
            }

            addGeneTask( $config, $def, $tasks, $target_dir, $choose_task, $choose_task, ".meta.rds", "cell_type", "seurat_cell_type" );

            my $pseudo_count_task = $seurat_task . "_multires" . get_next_index($def, $multiresKey) . "_pseudo_count";
            add_pseudo_count($config, $def, $tasks, $target_dir, $pseudo_count_task, $obj_ref, "seurat_cell_type");

            if(defined $def->{bubble_plots}){
              my $bubble_task = $seurat_task . "_multires" . get_next_index($def, $multiresKey) . "_bubblemap";
              add_bubble_plots($config, $def, $tasks, $target_dir, $bubble_task, $choose_task, undef, undef, undef, ".multi_choose_dot.html");
            }

            if ( $perform_comparison ) {
              if ( defined $def->{"DE_cluster_pairs"} ) {
                addEdgeRTask( $config, $def, $tasks, $target_dir, $choose_task, $choose_task, ".meta.rds", "cell_type", "seurat_cell_type", 1, 0, $DE_by_cell );
              }

              for my $deByOption (@deByOptions) {
                my $DE_by_celltype = $deByOption eq "DE_by_celltype";
                my $edgeRtask = addEdgeRTask( $config, $def, $tasks, $target_dir, $choose_task, $choose_task, ".meta.rds", "cell_type", "seurat_cell_type", 0, $DE_by_celltype, $DE_by_cell );
              }
            }
          }
        }
      }

      if(getValue($def, "perform_fix_resolution")){
        my $cluster_task = $seurat_task . "_cluster_res" . getValue( $def, "resolution" );
        $config->{$cluster_task} = {
          class                    => "CQS::UniqueR",
          perform                  => 1,
          target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $cluster_task,
          rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_cluster.r",
          parameterFile1_ref => [$seurat_task, ".rds"],
          parameterSampleFile1     => {
            Mtpattern             => getValue( $def, "Mtpattern" ),
            rRNApattern           => getValue( $def, "rRNApattern" ),
            Remove_rRNA        => getValue( $def, "Remove_rRNA" ),
            Remove_MtRNA        => getValue( $def, "Remove_MtRNA" ),
            resolution            => getValue( $def, "resolution" ),
            pca_dims              => getValue( $def, "pca_dims" ),
            by_sctransform        => getValue( $def, "by_sctransform" ),
            reduction             => $reduction,
          },
          output_file_ext      => ".final.rds",
          output_other_ext  => ".cluster.csv;.cluster.meanexp.csv;.umap.sample_cell.png;.cluster_sample.csv;.cluster_sample_percByCluster.png;.cluster_sample_percBySample.png",
          sh_direct            => 1,
          pbs                  => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "10",
            "mem"       => getValue($def, "seurat_mem")
          },
        };
        push( @$summary, $cluster_task );

        push(@report_files, ($cluster_task, ".umap.sample_cell.png", $cluster_task, ".cluster_sample.csv", $cluster_task, ".cluster_sample_percByCluster.png", $cluster_task, ".cluster_sample_percBySample.png" ));
        push(@report_names, "cluster_umap_sample_cell_png", "cluster_sample_cell_csv", "cluster_sample_percByCluster_png", "cluster_sample_percBySample_png");

        if(getValue($def, "perform_heatmap", 0)){
          my $heatmap_task = $cluster_task . "_heatmap";
          $config->{$heatmap_task} = {
            class                => "CQS::UniqueR",
            perform              => 1,
            target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $heatmap_task,
            rtemplate            => "../scRNA/scRNA_func.r,../scRNA/seurat_heatmap.r",
            parameterFile1_ref   => [ $cluster_task, ".final.rds" ],
            parameterFile2_ref   => [ $cluster_task, ".cluster.csv" ],
            parameterSampleFile1 => $def->{heatmap},
            output_file_ext      => ".heatmap.files",
            sh_direct            => 1,
            pbs                  => {
              "nodes"     => "1:ppn=1",
              "walltime"  => "10",
              "mem"       => getValue($def, "seurat_mem")
            },
          };
          push( @$summary, $heatmap_task );
        }

        $celltype_task = $cluster_task . "_celltype";
        if(getValue( $def, "annotate_tcell", 0)){
          getValue( $def, "HLA_panglao5_file");
          getValue( $def, "tcell_markers_file");
        }

        $config->{$celltype_task} = {
          class                    => "CQS::UniqueR",
          perform                  => 1,
          target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $celltype_task,
          rtemplate                => "../scRNA/scRNA_func.r,../scRNA/celltype_annotation.r",
          parameterFile1_ref =>  [$cluster_task, ".cluster.meanexp.csv"],
          parameterFile2_ref =>  [$cluster_task, ".cluster.csv"],
          parameterFile3_ref =>  [$cluster_task, ".final.rds"],
          parameterSampleFile1    => {
            species               => getValue( $def, "species" ),
            db_markers_file       => getValue( $def, "markers_file" ),
            curated_markers_file  => getValue( $def, "curated_markers_file", "" ),
            annotate_tcell        => getValue( $def, "annotate_tcell", 0),
            remove_subtype        => getValue( $def, "remove_subtype", ""),
            HLA_panglao5_file     => getValue( $def, "HLA_panglao5_file", "" ),
            tcell_markers_file    => getValue( $def, "tcell_markers_file", ""),
            by_sctransform        => getValue( $def, "by_sctransform" ),
            summary_layer_file => $def->{summary_layer_file},
          },
          output_file_ext      => ".celltype.csv",
          output_other_ext  => ".celltype_cluster.csv;.celltype.rds;.summary_layer.png;.cluster.png;.celltype.group.label.png;.score.png",
          sh_direct            => 1,
          pbs                  => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "40gb"
          },
        };

        if(defined $config->{groups}){
          $config->{$celltype_task}{parameterSampleFile2_ref} = "groups";
        }

        push( @$summary, $celltype_task );
        my $celltype_file     = ".celltype.csv";
        my $celltype_cluster_file      = ".celltype_cluster.csv";
        $celltype_name     = "cellactivity_clusters";
        my $cluster_name      = "seurat_cellactivity_clusters";

        push(@report_files, ($cluster_task, ".final.rds", 
                            $celltype_task, ".celltype.csv", 
                            $celltype_task, ".celltype.rds",
                            $celltype_task, ".summary_layer.png", 
                            $celltype_task, ".cluster.png",
                            $celltype_task, ".celltype.group.label.png",
                            $celltype_task, ".score.png"
                            ));
        push(@report_names, ("seurat_rds", "activity_celltype", "activity_rds", "activity_summary_layer_png", "activity_cluster_png", "activity_group_png", "activity_score_png"));

        my $find_markers = $cluster_task . "_celltype_markers";
        $config->{$find_markers} = {
          class                => "CQS::UniqueR",
          perform              => 1,
          target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $find_markers,
          rtemplate            => "../scRNA/scRNA_func.r,../scRNA/celltype_markers.r",
          parameterFile1_ref   => [ $cluster_task, ".final.rds" ],
          parameterFile2_ref   => [ $cluster_task, ".cluster.csv" ],
          parameterFile3_ref   => [ $celltype_task, ".celltype.csv" ],
          parameterSampleFile1 => {
            by_sctransform        => getValue( $def, "by_sctransform" ),
            celltype_name         => $celltype_name
          },
          output_file_ext      => "_celltype.all_display_markers.csv",
          output_other_ext     => "_celltype.all_max_markers.csv;_celltype.all_figures.csv;.top10.heatmap.png",
          sh_direct            => 1,
          pbs                  => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        };
        push( @$summary, $find_markers );

        push(@report_files, ($find_markers, "_celltype.all_display_markers.csv", $find_markers, "_celltype.all_figures.csv", $find_markers, ".top10.heatmap.png"));
        push(@report_names, ("celltype_markers_csv", "celltype_markers_png_list", "celltype_markers_heatmap_png"));

        if(getValue($def, "plot_gsea_genes", 0)){
          my $gesa_genes_task  = $celltype_task . "_gesa_genes";
          $config->{$gesa_genes_task} = {
            class              => "CQS::UniqueR",
            perform            => 1,
            target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $gesa_genes_task,
            rtemplate          => "../scRNA/scRNA_func.r;../scRNA/gsea_genes.r",
            parameterSampleFile1 => {
              by_sctransform        => getValue( $def, "by_sctransform" ),
            },
            parameterFile1_ref => [ $cluster_task, ".final.rds" ],
            parameterFile2_ref   => [ $celltype_task, ".celltype.csv" ],
            parameterFile3     => getValue($def, "gsea_genes_gmt"),
            parameterSampleFile2 => $def->{groups},
            output_file_ext    => ".figure.files.csv",
            sh_direct          => 1,
            pbs                => {
              "nodes"     => "1:ppn=1",
              "walltime"  => "1",
              "mem"       => "10gb"
            },
          };
          push( @$summary, $gesa_genes_task );
        }

        if (getValue( $def, "perform_scMRMA", 0 ) ) {
          addScMRMA( $config, $def, $tasks, $target_dir, $project_name, $cluster_task );
        }

        my $parameterSampleFile5_ref = undef;
        if (getValue( $def, "perform_CHETAH", 0 ) ) {
          my $CHETAH_name= $cluster_task . "_CHETAH";
          addCHETAH( $config, $def, $tasks, $target_dir, $project_name, $CHETAH_name, $cluster_task );
          push @report_files, ($CHETAH_name, ".CHETAH.png");
          push @report_names, "chetah_png";
        }

        addGeneTask( $config, $def, $tasks, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name );

        if ( getValue( $def, "perform_rename_cluster" ) ) {
          my $rename_cluster = $celltype_task . "_rename";
          $config->{$rename_cluster} = {
            class                => "CQS::UniqueR",
            perform              => 1,
            target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $rename_cluster,
            rtemplate            => "../scRNA/scRNA_func.r,../scRNA/rename_cluster.r",
            parameterFile1_ref   => [ $cluster_task, ".final.rds" ],
            parameterFile2_ref   => [ $celltype_task, ".cluster.csv" ],
            parameterFile3_ref   => [ $celltype_task, ".celltype.csv" ],
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

          $celltype_task = $rename_cluster;
          $celltype_file     = ".rename_celltype.csv";
          $celltype_cluster_file  = ".rename_cluster.csv";
          $celltype_name     = "renamed_cellactivity_clusters";
          $cluster_name      = "seurat_renamed_cellactivity_clusters";

          addGeneTask( $config, $def, $tasks, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name );
        }


        if(defined $def->{bubblemap_file}){
          my $bubblemap_name = $celltype_task . "_bubblemap";
          $config->{$bubblemap_name} = {
            class                => "CQS::UniqueR",
            perform              => 1,
            target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $bubblemap_name,
            rtemplate            => "../scRNA/scRNA_func.r,../scRNA/seurat_bubblemap.r",
            parameterFile1_ref   => [ $cluster_task, ".final.rds" ],
            parameterFile2_ref   => [ $cluster_task, ".cluster.csv" ],
            parameterFile3_ref   => [ $celltype_task, ".celltype.csv" ],
            parameterFile4       => $def->{bubblemap_file},
            parameterSampleFile1 => {
              by_sctransform        => getValue( $def, "by_sctransform" ),
              celltype_name => $celltype_name,
              cluster_name => $cluster_name,
              bubblemap_use_order => getValue($def, "bubblemap_use_order", 0),
            },
            output_file_ext      => ".bubblemap.png",
            sh_direct            => 1,
            pbs                  => {
              "nodes"     => "1:ppn=1",
              "walltime"  => "1",
              "mem"       => "10gb"
            },
          };
          push( @$summary, $bubblemap_name );

          push(@report_files, ($bubblemap_name, ".bubblemap.png"));
          push(@report_names, "bubblemap_png");
        }

        if ( getValue( $def, "perform_recluster" ) ) {
          my $recluster_task = $cluster_task . "_recluster";
          $config->{$recluster_task} = {
            class                => "CQS::UniqueR",
            perform              => 1,
            target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $recluster_task,
            rtemplate            => "../scRNA/scRNA_func.r;../scRNA/seurat_recluster.r",
            parameterFile1_ref   => [ $cluster_task, ".final.rds" ],
            parameterFile2_ref => [ $celltype_task, $celltype_cluster_file ],
            parameterFile3_ref => [ $celltype_task, ".rds" ],
            parameterSampleFile1 => $def->{recluster},
            parameterSampleFile2 => {
              recluster_celltypes => getValue( $def, "recluster_celltypes", "" ),
              Mtpattern           => getValue( $def, "Mtpattern" ),
              rRNApattern         => getValue( $def, "rRNApattern" ),
              Remove_rRNA         => getValue( $def, "Remove_rRNA" ),
              Remove_MtRNA        => getValue( $def, "Remove_MtRNA" ),
              resolution          => getValue( $def, "recluster_resolution" , getValue( $def, "resolution" )),
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
              "nodes"     => "1:ppn=1",
              "walltime"  => "1",
              "mem"       => "40gb"
            },
          };
          push( @$summary, $recluster_task );
        }

        if ( $def->{perform_antibody_vis} ) {
          addAntibody( $config, $def, $tasks, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name );
        }

        if(defined $clonotype_convert){
          addClonotypeVis($config, $def, $tasks, $target_dir, $clonotype_convert . "_vis", [ $cluster_task, ".final.rds" ], [ $celltype_task, $celltype_cluster_file ], $clonotype_convert);

          if(defined $celltype_task){
            my $clonetype_cluster = $clonotype_db;
            $clonetype_cluster =~ s/5_db/6_cluster/ig;
            $clonetype_cluster =~ s/6_db/7_cluster/ig;
            $config->{$clonetype_cluster} = {
              class                      => "CQS::UniqueR",
              perform                    => 1,
              target_dir                 => $target_dir . "/" . $clonetype_cluster,
              rtemplate                  => "../scRNA/clonotype_cluster.r",
              output_to_result_directory => 1,
              output_file_ext            => ".clonotype_cluster.csv",
              parameterFile1_ref         => [ $clonotype_db ],
              parameterFile2_ref         => [ $celltype_task, $celltype_cluster_file ],
              sh_direct                  => 1,
              pbs                        => {
                "nodes"     => "1:ppn=1",
                "walltime"  => "23",
                "mem"       => "10gb"
              },
            };
            push(@$summary, $clonetype_cluster);
          }
        }

        if ( $def->{perform_marker_dotplot} ) {
          my $biomarker_dotplot_task  = $cluster_task . "_biomarker_dotplot";
          my $marker_dotplot_clusters = getValue( $def, "marker_dotplot_clusters", {"all" => ["all"]});
          $config->{$biomarker_dotplot_task} = {
            class              => "CQS::UniqueR",
            perform            => 1,
            target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $biomarker_dotplot_task,
            rtemplate          => "../scRNA/scRNA_func.r;../scRNA/biomarker_dotplot.r",
            source             => $marker_dotplot_clusters,
            parameterFile1_ref => [ $cluster_task, ".final.rds" ],
            parameterFile2_ref => [ $celltype_task, $celltype_cluster_file ],
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
          my $curated_gene_dotplot_task  = $celltype_task . "_curated_gene_dotplot";

          my $curated_gene_dotplot = $def->{curated_gene_dotplot};
          my $curated_gene_files = $def->{curated_gene_files};

          if(not defined $curated_gene_dotplot and not defined $curated_gene_files) {
            die "Define curated_gene_dotplot or curated_gene_files first.";
          }

          if (not defined $curated_gene_dotplot){
            $curated_gene_dotplot = {};
          }

          if(defined $curated_gene_files){
            for my $cf (keys %$curated_gene_files){
              my $cf_file = $curated_gene_files->{$cf};
              my @cf_genes = read_file($cf_file, chomp => 1);
              my @pure_genes = grep { /\S/ } @cf_genes;
              $curated_gene_dotplot->{$cf} = {
                clusters => ["all"],
                genes => \@pure_genes
              };
            }
          }

          #print(Dumper($curated_gene_dotplot));

          my ($expanded_gene_def, $clusters, $genes) = parse_curated_genes($curated_gene_dotplot);

          $config->{$curated_gene_dotplot_task} = {
            class              => "CQS::UniqueR",
            perform            => 1,
            target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $curated_gene_dotplot_task,
            rtemplate          => "../scRNA/scRNA_func.r;../scRNA/curated_gene_dotplot.r",
            source             => $expanded_gene_def,
            parameterFile1_ref => [ $cluster_task, ".final.rds" ],
            parameterFile2_ref => [ $celltype_task, $celltype_cluster_file ],
            parameterSampleFile1 => $genes,
            parameterSampleFile2 => $clusters,
            parameterSampleFile3 => {
              cluster_name => "seurat_clusters",
              sort_cluster_name => $celltype_name,
              display_cluster_name => $cluster_name,
              by_sctransform => getValue($def, "by_sctransform"),
            },
            parameterSampleFile4_ref   => "groups",
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
          my $tcell_clusters_task  = $celltype_task . "_tcells";
          $config->{$tcell_clusters_task} = {
            class              => "CQS::UniqueR",
            perform            => 1,
            target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $tcell_clusters_task,
            rtemplate          => "../scRNA/extractTcells.r",
            parameterFile1_ref => [ $cluster_task, ".final.rds" ],
            parameterFile2     => $def->{marker_genes_file},
            parameterFile3_ref => [ $celltype_task, $celltype_cluster_file ],
            output_file_ext    => ".count.files.csv",
            rCode              => "celltype_name='" . $celltype_name . "'; cluster_name='" . join( ',', @$tcell_clusters ) . "'",
            sh_direct          => 1,
            pbs                => {
              "nodes"     => "1:ppn=1",
              "walltime"  => "1",
              "mem"       => "10gb"
            },
          };
          push( @$summary, $tcell_clusters_task );
        }

        if ( $perform_comparison ) {
          if ( defined $def->{"DE_cluster_pairs"} ) {
            addEdgeRTask( $config, $def, $tasks, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name, 1, 0, $DE_by_cell );
          }

          for my $deByOption (@deByOptions) {
            my $DE_by_celltype = $deByOption eq "DE_by_celltype";
            addEdgeRTask( $config, $def, $tasks, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name, 0, $DE_by_celltype, $DE_by_cell );
          }
        }
      }

      if(getValue($def, "perform_report", 1)){
        my $report_task= "report";
        my $additional_rmd_files = "Functions.Rmd;../scRNA/scRNA_func.r";
        my $report_options = merge({
              prefix => $project_name,
              summary_layer_file => $def->{summary_layer_file},
              celltype_name => $celltype_name
            }, merge($config->{$seurat_task}{parameterSampleFile1}, $config->{$celltype_task}{parameterSampleFile1}));
        #print(Dumper($report_options));

        $config->{$report_task} = {
          class                    => "CQS::BuildReport",
          perform                  => 1,
          target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $report_task,
          report_rmd_file          => "../scRNA/report.rmd",
          additional_rmd_files     => $additional_rmd_files,
          parameterSampleFile1_ref => \@report_files,
          parameterSampleFile1_names => \@report_names,
          parameterSampleFile2 => $report_options,
          parameterSampleFile3 => [],
          output_file_ext      => "_ur.html",
          sh_direct            => 1,
          pbs                  => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        };
        if(defined $hto_task){
          $config->{$report_task}{parameterSampleFile4_ref} = [$hto_task, ".HTO.class.dist.png"];
        }
        push( @$summary, $report_task );
      }
    }
  }

  if(defined $def->{seurat_object_file}){
    if ( $def->{perform_comparison} ) {
      if ( defined $def->{"DE_cluster_pairs"} ) {
        addEdgeRTask( $config, $def, $tasks, $target_dir, $def->{seurat_object_file}, undef, undef, getValue($def, "DE_cluster_name"), undef, 1, 0, 1 );
      }
    }
  }

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => $tasks,
    },
    sh_direct => 0,
    cluster   => $cluster,
    pbs       => {
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => $def->{sequencetask_run_time},
      "mem"       => "40gb"
    },
  };

  if($cluster ne "slurm"){
    $config->{sequencetask}{source}{step2} = $summary;
  }

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
