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
  initDefaultValue( $def, "hemoglobinPattern",   "^HB[^(P)]|^Hb[^(p)]" );
  
  initDefaultValue( $def, "Remove_rRNA",         1 );
  initDefaultValue( $def, "Remove_MtRNA",        1 );
  initDefaultValue( $def, "regress_by_percent_mt",   getValue($def, "Remove_MtRNA") ? 0 : 1 );
  initDefaultValue( $def, "Remove_hemoglobin",   0 );
  
  initDefaultValue( $def, "nFeature_cutoff_min", 200 );
  initDefaultValue( $def, "nFeature_cutoff_max", 10000 );
  initDefaultValue( $def, "nCount_cutoff",       500 );
  initDefaultValue( $def, "mt_cutoff",           20 );
  initDefaultValue( $def, "resolution",          0.5 );
  initDefaultValue( $def, "details_rmd",         "" );

  initDefaultValue( $def, "by_integration",        0 );
  initDefaultValue( $def, "by_sctransform",        0 );

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

  initDefaultValue( $def, "DE_by_sample_filter_minTPM",         1 );
  initDefaultValue( $def, "DE_by_sample_filter_cellPercentage", 0.5 );

  initDefaultValue( $def, "perform_webgestalt", 0 );

  initDefaultValue( $def, "perform_CHETAH", 0 );

  initDefaultValue( $def, "perform_dynamic_cluster", 1 );

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
  initDefaultValue( $def, "remove_subtype", "T cells,B cells,Fibroblasts,Neurons,Epithelial cells,Endothelial cells,Macrophages,Dendritic cells"),
  initDefaultValue( $def, "best_resolution_min_markers", 20);

  return $def;
}


sub getScRNASeqConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeScRNASeqDefaultOptions($def);

  my $project_name = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  if($cluster eq "slurm"){
    push(@$individual, @$summary);
    $summary = $individual;
  }


  my $target_dir      = $def->{target_dir};
  my $groups_ref      = defined $def->{groups} ? "groups" : undef;
  my $aligner         = $def->{aligner};
  my $star_option     = $def->{star_option};
  my $count_table_ref = "files";

  $config->{bam_files} = $def->{bam_files};

  my $perform_split_hto_samples = getValue($def, "perform_split_hto_samples", 0);
  if($perform_split_hto_samples){
    my $hto_samples = getValue($def, "HTO_samples");
    $def->{hto_file_names} = [sort keys %$hto_samples];
  }
  my $hto_file_names = $def->{hto_file_names};

  my $perform_clonotype_analysis = getValue($def, "perform_clonotype_analysis", 0);

  my $clonotype_4_convert = undef;
  if ($perform_clonotype_analysis){
    if ((not defined $def->{files}) || (not $perform_split_hto_samples)) {
      $config->{vdj_json_files} = getValue($def, "vdj_json_files");
      addClonotypeMerge($config, $def, $summary, $target_dir, "clonotype_1_merge", ["vdj_json_files", "all_contig_annotations.json"]);
      addEnclone($config, $def, $summary, "clonotype_2_enclone", $target_dir, ["clonotype_1_merge", ".json\$"] );
      $clonotype_4_convert = addEncloneToClonotype($config, $def, $summary, $target_dir, "clonotype_3_convert", "clonotype_2_enclone", ["clonotype_1_merge", ".cdr3\$"]);
    }
  }

  my $hla_merge = undef;
  my $individual_qc_task = "individual_qc";
  my $qc_filter_config_file = $target_dir . "/qc_filter_config.txt";
  if (defined $def->{files}){
    my @report_files = ();
    my @report_names = ();
    my $hto_task = undef;
    my $hto_ref = undef;
    my $hto_sample_file = undef;
    my $hto_summary_task = undef;

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

    if( $perform_split_hto_samples ) {
      my $preparation_task = "hto_samples_preparation";
      $config->{$preparation_task} = {
        class => "CQS::UniqueR",
        target_dir => "${target_dir}/$preparation_task",
        rtemplate => "../scRNA/split_samples_utils.r,../scRNA/split_samples_preparation.r",
        option => "",
        parameterSampleFile1_ref => $hto_file_ref,
        parameterSampleFile2 => {
          hto_regex => getValue($def, "hto_regex", ""),
          nFeature_cutoff_min => getValue($def, "nFeature_cutoff_min"),
          hto_non_zero_percentage => getValue($def, "hto_non_zero_percentage", 0.2),
        },
        output_perSample_file => "parameterSampleFile1",
        output_perSample_file_byName => 1,
        output_perSample_file_ext => ".hto.rds;.barcodes.tsv",
        output_to_same_folder => 1,
        can_result_be_empty_file => 0,
        sh_direct   => 1,
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary, $preparation_task );
      $hto_file_ref = [ $preparation_task, ".hto.rds"];

      my $r_script = undef;
      if ( getValue($def, "split_hto_samples_by_cutoff", 0) ) {
        if(getValue($def, "use_cutoff_v2", 0)){
          $r_script = "../scRNA/split_samples_utils.r,../scRNA/split_samples_cutoff_all_v2.r";
          $hto_task = "hto_samples_cutoff_all_v2";
        }else{
          $r_script = "../scRNA/split_samples_utils.r,../scRNA/split_samples_cutoff_all.r";
          $hto_task = "hto_samples_cutoff_all";
        }
      } else {
        $r_script = "../scRNA/split_samples_utils.r,../scRNA/split_samples_seurat_all.r";
        $hto_task = "hto_samples_HTODemux_all";
      }

      $config->{$hto_task} = {
        class => "CQS::UniqueR",
        target_dir => "${target_dir}/$hto_task",
        rtemplate => $r_script,
        option => "",
        parameterSampleFile1_ref => $hto_file_ref,
        #parameterSampleFile2 => $def->{split_hto_samples_cutoff_point},
        parameterSampleFile3 => {
          hto_ignore_exists => getValue($def, "hto_ignore_exists", 0),
          cutoff_file => getValue($def, "cutoff_file", ""),
          umap_min_dist => getValue($def, "hto_umap_min_dist", 0.3),
          umap_num_neighbors => getValue($def, "hto_umap_num_neighbors", 30),
        },
        output_perSample_file => "parameterSampleFile1",
        output_perSample_file_byName => 1,
        output_perSample_file_ext => ".HTO.umap.class.png;.HTO.csv;.HTO.data.csv;.HTO.umap.rds",
        sh_direct   => 1,
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary, $hto_task );
      $hto_ref = [ $hto_task, ".HTO.csv" ];

      $hto_summary_task = $hto_task . "_summary";
      $config->{$hto_summary_task} = {
        class => "CQS::UniqueR",
        target_dir => "${target_dir}/${hto_summary_task}",
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
      push( @$summary, $hto_summary_task );

      push (@report_files, ($hto_summary_task, ".HTO.summary.global.png"));
      push (@report_names, "hto_summary_png");

      my $perform_souporcell = getValue($def, "perform_souporcell", 0);
      if($perform_souporcell) {
        my $fasta = getValue($def, "fasta_file");
        my $skip_remap = getValue($def, "skip_remap", 0);
        my $common_variants = getValue($def, "common_variants", "");

        my $hto_souporcell_task = "hto_souporcell" . ($skip_remap ? "_skip_remap" : "_remap");

        my $skip_remap_option = $skip_remap ? "--skip_remap SKIP_REMAP" : "";
        my $common_variants_option = ($common_variants eq "") ? "" : "--common_variants $common_variants";
        my $souporcell_thread = getValue($def, "souporcell_cpu", "16");

        $config->{$hto_souporcell_task} = {
          class => "CQS::ProgramWrapperOneToOne",
          target_dir => "${target_dir}/$hto_souporcell_task",
          interpretor => "",
          program => "souporcell_pipeline.py",
          check_program => 0,
          docker_prefix => "souporcell_",
          option => "-i __FILE__ -b __FILE2__ -f $fasta -t $souporcell_thread -o . -k __FILE3__ $common_variants_option $skip_remap_option
          
#__OUTPUT__
",
          source_arg => "-i",
          source_ref => "bam_files",
          parameterSampleFile2_arg => "-b",
          parameterSampleFile2_ref => [ $preparation_task, ".barcodes.tsv"],
          parameterSampleFile3_arg => "-k",
          parameterSampleFile3 => getValue($def, "souporcell_tag_number"),
          output_arg => "-o",
          output_file_prefix => "",
          output_no_name => 1,
          output_file_ext => "clusters.tsv",
          output_to_same_folder => 0,
          can_result_be_empty_file => 0,
          use_tmp_folder => getValue($def, "use_tmp_folder", 1),
          sh_direct   => 0,
          pbs => {
            "nodes"     => "1:ppn=" . $souporcell_thread,
            "walltime"  => getValue($def, "souporcell_walltime", "47"),
            "mem"       => getValue($def, "souporcell_mem", "40gb")
          },
        };
        push( @$individual, $hto_souporcell_task );

        my $hto_integration_task = $hto_souporcell_task . "_cutoff_integration";
        $config->{$hto_integration_task} = {
          class => "CQS::UniqueR",
          target_dir => "${target_dir}/${hto_integration_task}",
          rtemplate => "../scRNA/hto_soupercell_integration.r",
          option => "",
          parameterSampleFile1_ref => $hto_souporcell_task,
          parameterSampleFile2_ref => $hto_ref,
          parameterSampleFile3_ref => [ $hto_task, ".umap.rds" ],
          parameterSampleFile4     => $def->{ignore_souporcell_cluster},
          output_perSample_file => "parameterSampleFile1",
          output_perSample_file_byName => 1,
          output_perSample_file_ext => ".HTO.png;.meta.rds",
          sh_direct   => 1,
          pbs => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        };
        push( @$summary, $hto_integration_task );

        $hto_ref = [ $hto_integration_task, ".meta.rds" ];
      }

      if(defined $def->{HTO_samples}){
        $hto_sample_file = write_HTO_sample_file($def);
      }

      my $perform_arcasHLA = getValue($def, "perform_arcasHLA", 0);

      if($perform_arcasHLA){
        if ( not defined $def->{bam_files}){
          die "Define bam_files for perform_arcasHLA";
        }

        if (not defined $def->{HTO_samples}) {
          die "Define HTO_samples for split bam files";
        }

        $config->{HTO_samples} = $def->{HTO_samples};
        $config->{bam_files} = $def->{bam_files};

        my $hto_bam_task = "hto_bam";
        $config->{$hto_bam_task} = {
          class => "CQS::ProgramWrapperOneToManyFile",
          target_dir => "${target_dir}/hto_bam",
          interpretor => "python3",
          program => "../scRNA/split_samples.py",
          check_program => 1,
          option => "-o .",
          source_arg => "-s",
          source_ref => ["HTO_samples"],
          parameterSampleFile2_arg => "-b",
          parameterSampleFile2_ref => ["bam_files"],
          parameterSampleFile3_arg => "-i",
          parameterSampleFile3_ref => $hto_ref,
          output_arg => "-o",
          output_file_prefix => "",
          output_file_key => 0,
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
        push( @$individual, $hto_bam_task );

        $hla_merge = addArcasHLA($config, $def, $individual, $target_dir, $project_name, $hto_bam_task . "_", $hto_bam_task);        
      }

      if($perform_clonotype_analysis){
        if (not defined $def->{HTO_samples}) {
          die "Define HTO_samples for split vdj json files";
        }

        $config->{HTO_samples} = $def->{HTO_samples};
        $config->{vdj_json_files} = getValue($def, "vdj_json_files");
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
      if(getValue($def, "perform_arcasHLA", 0)){
        getValue($def, "bam_files");
        if (defined $def->{singleend_bam_files}){
          $config->{singleend_bam_files} = $def->{singleend_bam_files};
          $hla_merge = addArcasHLA($config, $def, $individual, $target_dir, $project_name, "", "bam_files", "singleend_bam_files");        
        }else{
          $hla_merge = addArcasHLA($config, $def, $individual, $target_dir, $project_name, "", "bam_files");        
        }
      }
    }

    my $clonotype_db = undef;

    #die "clonotype_4_convert=$clonotype_4_convert";
    #die "hla_merge=$hla_merge";

    if(defined $clonotype_4_convert) {
      $clonotype_db = $clonotype_4_convert . "_db";
      addClonotypeDB($config, $def, $summary, $target_dir, $clonotype_db, $clonotype_4_convert);
    }

    if($def->{perform_individual_qc}){
      my $qc_files_ref;
      if(defined $def->{qc_files}){
        $qc_files_ref = "qc_files";
        $config->{qc_files} = $def->{qc_files};
      }else{
        $qc_files_ref = "files";
      }
      $config->{$individual_qc_task} = {
        class => "CQS::UniqueRmd",
        target_dir => "${target_dir}/$individual_qc_task",
        report_rmd_file => "../scRNA/individual_qc.Rmd",
        additional_rmd_files => "../scRNA/markerCode_filter.R;../scRNA/scRNA_func.r",
        option => "",
        parameterSampleFile1_ref => $qc_files_ref,
        parameterSampleFile2 => {
          species => getValue($def, "species"),
          Mtpattern             => getValue( $def, "Mtpattern" ),
          rRNApattern           => getValue( $def, "rRNApattern" ),
          Remove_rRNA        => getValue( $def, "Remove_rRNA" ),
          Remove_MtRNA        => getValue( $def, "Remove_MtRNA" ),
          nFeature_cutoff_min   => getValue( $def, "nFeature_cutoff_min" ),
          nFeature_cutoff_max   => getValue( $def, "nFeature_cutoff_max" ),
          nCount_cutoff         => getValue( $def, "nCount_cutoff" ),
          mt_cutoff             => getValue( $def, "mt_cutoff" ),
          species               => getValue( $def, "species" ),
          resolution            => getValue( $def, "resolution" ),
          pca_dims              => getValue( $def, "pca_dims" ),
          markers_file     => getValue( $def, "markers_file" ),
          curated_markers_file     => getValue( $def, "curated_markers_file", "" ),
          bubblemap_file => getValue( $def, "bubblemap_file", "" ),
        },
        output_file_ext => "objectlist.Rdata",
        output_no_name => 1,
        can_result_be_empty_file => 0,
        sh_direct   => 1,
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };

      if($perform_split_hto_samples){
        $config->{$individual_qc_task}{parameterSampleFile3_ref} = $hto_ref;
        $config->{$individual_qc_task}{parameterSampleFile2}{hto_sample_file} = $hto_sample_file;
      }

      if( ! -e $qc_filter_config_file){
        open(my $qc, '>', $qc_filter_config_file) or die $!;
        print $qc "sample\tnFeature_cutoff_min\tnFeature_cutoff_max\tnCount_cutoff\tmt_cutoff\tcluster_remove\n";
        my $files = $def->{files};
        for my $fname (sort keys %$files){
          print $qc "$fname\t" . 
                    getValue( $def, "nFeature_cutoff_min" ) . "\t" . 
                    getValue( $def, "nFeature_cutoff_max" ) . "\t" . 
                    getValue( $def, "nCount_cutoff" ) . "\t" .
                    getValue( $def, "mt_cutoff" ) . "\t\n";
        }
        close($qc);
      }
      
      push( @$summary, $individual_qc_task );
    }

    if ( getValue( $def, "perform_scRNABatchQC" ) ) {
      $config->{scRNABatchQC} = {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $target_dir . "/" . getNextFolderIndex($def) . "scRNABatchQC",
        rtemplate                => "../scRNA/scRNABatchQC.r",
        parameterSampleFile1_ref => "files",
        rCode                    => "webgestalt_organism='" . getValue( $def, "webgestalt_organism" ) . "'",
        output_file_ext      => ".html",
        sh_direct                => 1,
        pbs                      => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary, "scRNABatchQC" );
    }

    my $essential_gene_task = "essential_genes";
    $config->{$essential_gene_task} = {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $essential_gene_task,
      rtemplate                => "../scRNA/scRNA_func.r,../scRNA/essential_genes.r",
      parameterSampleFile1     => {
        species               => getValue( $def, "species" ),
        db_markers_file       => getValue( $def, "markers_file" ),
        curated_markers_file  => getValue( $def, "curated_markers_file", "" ),
        HLA_panglao5_file     => getValue( $def, "HLA_panglao5_file", "" ),
        remove_subtype        => getValue( $def, "remove_subtype", ""),
        bubblemap_file        => $def->{bubblemap_file},
      },
      parameterSampleFile2    => get_marker_gene_dict($def),
      rCode                    => "",
      output_file_ext          => ".txt",
      sh_direct                => 1,
      pbs                      => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    };
    push( @$summary, $essential_gene_task );

    my $seurat_task;
    if ( getValue( $def, "perform_seurat" ) ) {
      my $reduction = "pca";

      if (getValue($def, "perform_seurat_oldversion", 0)){
        $seurat_task = "seurat" . ( getValue( $def, "by_sctransform" ) ? "_sct" : "" ) . ( getValue( $def, "by_integration" ) ? "_igr" : "" ) . ( getValue( $def, "pool_sample" ) ? "_pool" : "" );
        $config->{$seurat_task} = {
          class                    => "CQS::UniqueR",
          perform                  => 1,
          target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $seurat_task,
          rtemplate                => "../scRNA/scRNA_func.r,../scRNA/analysis.rmd",
          parameterSampleFile1_ref => "files",
          parameterSampleFile2     => {
            Mtpattern             => getValue( $def, "Mtpattern" ),
            rRNApattern           => getValue( $def, "rRNApattern" ),
            Remove_rRNA        => getValue( $def, "Remove_rRNA" ),
            Remove_MtRNA        => getValue( $def, "Remove_MtRNA" ),
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
            integration_by_harmony=> getValue( $def, "integration_by_harmony", 1),
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
            "walltime"  => "12",
            "mem"       => "40gb"
          },
        };
        push( @$summary, $seurat_task );
      }
      
      my $seurat_rawdata = "seurat_rawdata";
      $config->{$seurat_rawdata} = {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $seurat_rawdata,
        rtemplate                => "../scRNA/scRNA_func.r;../scRNA/seurat_rawdata.r",
        parameterSampleFile1_ref => "files",
        parameterSampleFile2     => {
          Mtpattern             => getValue( $def, "Mtpattern" ),
          rRNApattern           => getValue( $def, "rRNApattern" ),
          hemoglobinPattern     => getValue( $def, "hemoglobinPattern" ),
          species               => getValue( $def, "species" ),
          pool_sample           => getValue( $def, "pool_sample" ),
          hto_sample_file       => $hto_sample_file,
        },
        parameterSampleFile3 => $def->{"pool_sample_groups"},
        parameterSampleFile4_ref => $hto_ref,
        output_file_ext      => ".rawobj.rds",
        sh_direct            => 1,
        pbs                  => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "12",
          "mem"       => "40gb"
        },
      };
      if($def->{perform_individual_qc}){
        $config->{$seurat_rawdata}{parameterFile1_ref} = [$individual_qc_task, "objectlist.Rdata"];
        $config->{$seurat_rawdata}{parameterFile2} = $qc_filter_config_file;
      }
      push( @$summary, $seurat_rawdata );

      push (@report_files, ($seurat_rawdata, "rawobj.rds"));
      push (@report_names, "raw_obj");

      my @sample_names = keys %{$def->{files}};
      my $nsamples = scalar(@sample_names);
      my $by_integration = $nsamples > 1 ? getValue( $def, "by_integration" ) : 0;
      my $sct_str = getValue( $def, "by_sctransform" ) ? "_sct":"";

      my $preprocessing_rscript;
      if($by_integration){
        my $integration_by_harmony = getValue( $def, "integration_by_harmony", 1);
        if($integration_by_harmony){
          $seurat_task = "seurat${sct_str}_harmony";
          $preprocessing_rscript = "../scRNA/seurat_harmony.r";
          $reduction = "harmony";
        }else{
          $seurat_task = "seurat${sct_str}_integration";
          $preprocessing_rscript = "../scRNA/seurat_integration.r";
        }
      }else{
        $seurat_task = "seurat${sct_str}_merge";
        $preprocessing_rscript = "../scRNA/seurat_merge.r";
      }

      $config->{$seurat_task} = {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $seurat_task,
        rtemplate                => "../scRNA/scRNA_func.r,$preprocessing_rscript",
        parameterFile1_ref => [$seurat_rawdata, ".rawobj.rds"],
        parameterFile2_ref => $essential_gene_task,
        parameterSampleFile1     => {
          Mtpattern             => getValue( $def, "Mtpattern" ),
          rRNApattern           => getValue( $def, "rRNApattern" ),
          Remove_rRNA           => getValue( $def, "Remove_rRNA" ),
          Remove_MtRNA          => getValue( $def, "Remove_MtRNA" ),
          regress_by_percent_mt => getValue( $def, "regress_by_percent_mt" ),
          nFeature_cutoff_min   => getValue( $def, "nFeature_cutoff_min" ),
          nFeature_cutoff_max   => getValue( $def, "nFeature_cutoff_max" ),
          nCount_cutoff         => getValue( $def, "nCount_cutoff" ),
          mt_cutoff             => getValue( $def, "mt_cutoff" ),
          species               => getValue( $def, "species" ),
          resolution            => getValue( $def, "resolution" ),
          pca_dims              => getValue( $def, "pca_dims" ),
          by_integration        => $by_integration,
          by_sctransform        => getValue( $def, "by_sctransform" ),
          batch_for_integration => getValue( $def, "batch_for_integration" ),
        },
        output_file_ext      => ".final.rds,.qc.1.png,.qc.2.png,.qc.3.png,.qc.4.png,.sample_cell.csv,.final.png",
        sh_direct            => 1,
        pbs                  => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "12",
          "mem"       => "40gb"
        },
      };
      if ($def->{batch_for_integration}){
        $config->{$seurat_task}{parameterSampleFile2} = getValue($def, "batch_for_integration_groups");
      }

      push( @$summary, $seurat_task );

      push (@report_files, ($seurat_task, ".final.png", 
        $seurat_task, ".qc.1.png", 
        $seurat_task, ".qc.2.png", 
        $seurat_task, ".qc.3.png", 
        $seurat_task, ".qc.4.png", 
        $seurat_task, ".sample_cell.csv"));
      push (@report_names, ("seurat_merge_png", "seurat_qc_1_png", "seurat_qc_2_png", "seurat_qc_3_png", "seurat_qc_4_png", 
        "sample_cell_csv"));

      my $localization_ref = [ $seurat_task, ".final.rds" ];

      my $df_task = $seurat_task . "_doublet_finder";
      if(getValue($def, "perform_doublet_finder", 0)){
        addDoubletFinder($config, $def, $summary, $target_dir, $df_task, [$seurat_task, ".rds"], undef );
      }

      my $signacX_name = $seurat_task . "_SignacX";
      if (getValue( $def, "perform_SignacX", 0 ) ) {
        addSignac_only( $config, $def, $summary, $target_dir, $project_name, $signacX_name, $seurat_task, $reduction );

        if(getValue($def, "perform_dynamic_cluster_signacX")){
          my $dynamicKey = "dynamic_cluster_singnacX";
          my $scDynamic_task = $seurat_task . get_next_index($def, $dynamicKey) . "_dynamic_signacX";

          addDynamicClusterSignacX($config, $def, $summary, $target_dir, $scDynamic_task, $seurat_task, $essential_gene_task, $reduction, $signacX_name);
        }
      }

      my $celltype_name = undef;
      my $celltype_task = undef;

      if(getValue($def, "perform_dynamic_cluster")){
        my $dynamicKey = "dynamic_cluster";
        my $scDynamic_task = $seurat_task . get_next_index($def, $dynamicKey) . "_dynamic";

        addDynamicCluster($config, $def, $summary, $target_dir, $scDynamic_task, $seurat_task, $essential_gene_task, $reduction);

        my $subcluster_task;
        my $can_check_doublet=0;
        if(getValue($def, "perform_subcluster_v2", 1)){
          $subcluster_task = $seurat_task . get_next_index($def, $dynamicKey) . "_subcluster.v2";

          addSubClusterV2($config, $def, $summary, $target_dir, $subcluster_task, $seurat_task, $scDynamic_task, $essential_gene_task, $reduction, $signacX_name);

          if (defined $def->{"celltype_subclusters_table"}){
            my $choose_task = $seurat_task . get_next_index($def, $dynamicKey) . "_choose_res";
            $config->{$choose_task} = {
              class                    => "CQS::UniqueR",
              perform                  => 1,
              target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $choose_task,
              rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_celltype_subcluster.v2.choose.r",
              parameterFile1_ref => [$seurat_task, ".rds"],
              parameterFile2_ref => [$scDynamic_task, ".meta.rds"],
              parameterFile3_ref => $essential_gene_task,
              parameterFile4_ref => [$subcluster_task, ".files.csv"],
              parameterSampleFile1     => {
                pca_dims              => getValue( $def, "pca_dims" ),
                by_sctransform        => getValue( $def, "by_sctransform" ),
                regress_by_percent_mt => getValue( $def, "regress_by_percent_mt" ),
                reduction             => $reduction,
                species               => getValue( $def, "species" ),
                db_markers_file       => getValue( $def, "markers_file" ),
                curated_markers_file  => getValue( $def, "curated_markers_file", "" ),
                annotate_tcell        => getValue( $def, "annotate_tcell", 0),
                #remove_subtype        => getValue( $def, "remove_subtype", ""),
                remove_subtype        => "",
                HLA_panglao5_file     => getValue( $def, "HLA_panglao5_file", "" ),
                tcell_markers_file    => getValue( $def, "tcell_markers_file", ""),
                bubblemap_file        => $def->{bubblemap_file},
                bubblemap_use_order   => getValue($def, "bubblemap_use_order", 0),
                summary_layer_file    => $def->{summary_layer_file},
                celltype_layer        => "layer4",
                output_layer          => "cell_type",
              },
              parameterSampleFile2 => $def->{"subcluster_ignore_gene_files"},
              parameterSampleFile3 => $def->{"celltype_subclusters_table"},
              output_file_ext      => ".meta.rds",
              output_other_ext  => ".final.rds,.umap.png",
              sh_direct            => 1,
              pbs                  => {
                "nodes"     => "1:ppn=1",
                "walltime"  => "12",
                "mem"       => "40gb"
              },
            };
            push( @$summary, $choose_task );
            $subcluster_task = $choose_task;
            $can_check_doublet=1;

            if(defined $clonotype_4_convert) {
              if( getValue($def, "perform_gliph2", 0) ) {
                my $prepare_task = (defined $hla_merge) ? "tcr_hla_data" : "tcr_data";
                #die "gliph2_task=$gliph2_task";
                $config->{$prepare_task} = {
                  class                => "CQS::UniqueR",
                  perform              => 1,
                  target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $prepare_task,
                  rtemplate            => "../scRNA/gliph2_prepare.r",
                  parameterFile1_ref   => [ $choose_task, ".meta.rds"],
                  parameterFile2_ref   => $clonotype_4_convert,
                  parameterSampleFile1 => getValue($def, "gliph2_config"),
                  parameterSampleFile2 => {
                    gliph2_hla_condition => getValue($def, "gliph2_hla_condition"),
                  },
                  output_file_ext      => ".tcr.CD4.txt,.tcr.CD8.txt",
                  output_other_ext     => "",
                  sh_direct            => 1,
                  pbs                  => {
                    "nodes"     => "1:ppn=1",
                    "walltime"  => "12",
                    "mem"       => "40gb"
                  },
                };
                if( defined $hla_merge ) {
                  $config->{$prepare_task}{parameterFile3_ref} = $hla_merge; 
                  $config->{$prepare_task}{output_other_ext} = ".hla.txt"; 
                }
                push( @$summary, $prepare_task );

                my $gliph2_config_file = dirname(__FILE__) . "/../scRNA/gliph2_config.txt";
                #die($gliph2_config_file);

                my $gliph2_task = $prepare_task . "_gliph2";
                $config->{$gliph2_task} = {
                  class                => "CQS::UniqueR",
                  perform              => 1,
                  target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $gliph2_task,
                  rtemplate            => "../scRNA/gliph2.r",
                  parameterFile1       => getValue($def, "gliph2_config_file", $gliph2_config_file),
                  parameterFile2_ref   => [ $prepare_task, '.tcr.CD4.txt'],
                  parameterFile3_ref   => [ $prepare_task, '.tcr.CD8.txt'],
                  parameterSampleFile1 => getValue($def, "gliph2_reference"),
                  output_file_ext      => "_HLA.txt",
                  sh_direct            => 1,
                  pbs                  => {
                    "nodes"     => "1:ppn=1",
                    "walltime"  => "12",
                    "mem"       => "40gb"
                  },
                };
                if( defined $hla_merge ) {
                  $config->{$gliph2_task}{parameterFile4_ref} = [$prepare_task, '.hla.txt'], 
                }
                push( @$summary, $gliph2_task );
              }
            }

            addComparison($config, $def, $summary, $target_dir, $choose_task, $choose_task, "", "cell_type", "seurat_cell_type");

            $localization_ref = [ $choose_task, ".final.rds" ];

            if(defined $def->{groups}){
              my $group_umap_task = $seurat_task . get_next_index($def, $dynamicKey) . "_group_umap";
              $config->{$group_umap_task} = {
                class                    => "CQS::UniqueR",
                perform                  => 1,
                target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $group_umap_task,
                rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_group_umap.r",
                parameterFile1_ref => [$choose_task, ".final.rds"],
                parameterSampleFile1     => $def->{groups},
                output_file_ext      => ".all.label.umap.png",
                sh_direct            => 1,
                pbs                  => {
                  "nodes"     => "1:ppn=1",
                  "walltime"  => "12",
                  "mem"       => "40gb"
                },
              };
              push( @$summary, $group_umap_task );
            }

            if(defined $clonotype_4_convert){
              addClonotypeVis($config, $def, $summary, $target_dir, $clonotype_4_convert . "_vis", [$choose_task, ".final.rds"], undef, $clonotype_4_convert);
            }

            if(defined $clonotype_db){
              addClonotypeCluster($config, $def, $summary, $target_dir, $clonotype_db . "_cluster", $clonotype_db, [$choose_task, ".meta.rds"]);
            }
          }
        }else{
          $subcluster_task = $seurat_task . get_next_index($def, $dynamicKey) . "_subcluster.v1";
          $config->{$subcluster_task} = {
            class                    => "CQS::UniqueR",
            perform                  => 1,
            target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $subcluster_task,
            rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_celltype_subcluster.r",
            parameterFile1_ref => [$seurat_task, ".rds"],
            parameterFile2_ref => [$scDynamic_task, ".meta.rds"],
            parameterFile3_ref => $essential_gene_task,
            parameterSampleFile1     => {
              pca_dims              => getValue( $def, "pca_dims" ),
              by_sctransform        => getValue( $def, "by_sctransform" ),
              regress_by_percent_mt => getValue( $def, "regress_by_percent_mt" ),
              reduction             => $reduction,
              species               => getValue( $def, "species" ),
              db_markers_file       => getValue( $def, "markers_file" ),
              curated_markers_file  => getValue( $def, "curated_markers_file", "" ),
              annotate_tcell        => getValue( $def, "annotate_tcell", 0),
              #remove_subtype        => getValue( $def, "remove_subtype", ""),
              remove_subtype        => "",
              HLA_panglao5_file     => getValue( $def, "HLA_panglao5_file", "" ),
              tcell_markers_file    => getValue( $def, "tcell_markers_file", ""),
              bubblemap_file        => $def->{bubblemap_file},
              bubblemap_use_order   => getValue($def, "bubblemap_use_order", 0),
              summary_layer_file    => $def->{summary_layer_file},
              celltype_layer        => "layer4",
              output_layer          => "cell_type",
              best_resolution_min_markers => getValue( $def, "best_resolution_min_markers" ),
            },
            parameterSampleFile2 => $def->{"subcluster_ignore_gene_files"},
            output_file_ext      => ".meta.rds",
            output_other_ext  => ".umap.png",
            sh_direct            => 1,
            pbs                  => {
              "nodes"     => "1:ppn=1",
              "walltime"  => "12",
              "mem"       => "40gb"
            },
          };
          push( @$summary, $subcluster_task );
          $can_check_doublet=1;
        }

        if($can_check_doublet){
          my $doublet_check_task = $seurat_task . get_next_index($def, $dynamicKey) . "_doublet";
          $config->{$doublet_check_task} = {
            class                    => "CQS::UniqueR",
            perform                  => 1,
            target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $doublet_check_task,
            rtemplate                => "countTableVisFunctions.R,../scRNA/scRNA_func.r,../scRNA/seurat_doublet_check.r",
            parameterFile1_ref       => [ $seurat_task, ".rds" ],
            parameterFile2_ref       => [ $subcluster_task, ".meta.rds" ],
            parameterFile3_ref       => [ $df_task, ".meta.rds" ],
            parameterFile4_ref       => [ $df_task, ".options.csv" ],
            parameterSampleFile1     => {
              cluster_layer           => "seurat_clusters",
              celltype_layer          => "cell_type",
              cluster_celltype_layer  => "seurat_cell_type",
              bubblemap_file        => $def->{bubblemap_file},
              by_sctransform        => getValue( $def, "by_sctransform" ),
            },
            output_file_ext      => ".doublet_perc.png",
            output_other_ext  => "",
            sh_direct            => 1,
            pbs                  => {
              "nodes"     => "1:ppn=1",
              "walltime"  => "12",
              "mem"       => "40gb"
            },
          };
          push( @$summary, $doublet_check_task );
        }

        if(getValue($def, "perform_dynamic_discard_cluster", 0)){
          my $discard_task = $seurat_task . get_next_index($def, $dynamicKey) . "_discard";
          $config->{$discard_task} = {
            class                    => "CQS::UniqueR",
            perform                  => 1,
            target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $discard_task,
            rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_discard_cluster.r",
            parameterFile1_ref => [$seurat_task, ".rds"],
            parameterFile2_ref => [$subcluster_task, ".meta.rds"],
            parameterSampleFile1     => {
              dynamic_discard_cluster => getValue( $def, "dynamic_discard_cluster" ),
              by_sctransform        => getValue( $def, "by_sctransform" ),
              reduction             => $reduction,
              bubblemap_file        => $def->{bubblemap_file},
              bubblemap_use_order   => getValue($def, "bubblemap_use_order", 0),
              cluster_layer           => "seurat_clusters",
              celltype_layer          => "cell_type",
              cluster_celltype_layer  => "seurat_cell_type",
            },
            output_file_ext      => ".meta.rds",
            output_other_ext  => ".umap.png,.meta.csv,.options.csv",
            sh_direct            => 1,
            pbs                  => {
              "nodes"     => "1:ppn=1",
              "walltime"  => "12",
              "mem"       => "40gb"
            },
          };
          push( @$summary, $discard_task );
          $subcluster_task = $discard_task;
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
        my $multires_task = $seurat_task . "_multires";
        $config->{$multires_task} = {
          class                    => "CQS::UniqueR",
          perform                  => 1,
          target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $multires_task,
          rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_multires.r",
          parameterFile1_ref => [$seurat_task, ".rds"],
          parameterFile2_ref => [$seurat_task, ".pca.txt"],
          parameterSampleFile1    => {
            Mtpattern             => getValue( $def, "Mtpattern" ),
            rRNApattern           => getValue( $def, "rRNApattern" ),
            Remove_rRNA           => getValue( $def, "Remove_rRNA" ),
            Remove_MtRNA          => getValue( $def, "Remove_MtRNA" ),
            pca_dims              => getValue( $def, "pca_dims" ),
            by_sctransform        => getValue( $def, "by_sctransform" ),
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
            plot_width            => getValue($def, "multires_plot_width", 6600),
            plot_height           => getValue($def, "multires_plot_height", 6000),
          },
          output_file_ext      => ".multires.rds",
          output_other_ext  => ".resolutions.csv;.umap.sample_cell.png",
          sh_direct            => 1,
          pbs                  => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "40gb"
          },
        };
        push( @$summary, $multires_task );

        if (getValue($def, "perform_multires_subcluster", 1)){
          my $recluster_task = $seurat_task . "_multires_subcluster";
          $config->{$recluster_task} = {
            class                    => "CQS::UniqueR",
            perform                  => 1,
            target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $recluster_task,
            rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_multires_subcluster.r",
            parameterFile1_ref => [$multires_task, ".rds"],
            parameterFile2_ref => [$multires_task, ".resolutions.csv"],
            parameterSampleFile1     => {
              Mtpattern             => getValue( $def, "Mtpattern" ),
              rRNApattern           => getValue( $def, "rRNApattern" ),
              Remove_rRNA           => getValue( $def, "Remove_rRNA" ),
              Remove_MtRNA          => getValue( $def, "Remove_MtRNA" ),
              pca_dims              => getValue( $def, "pca_dims" ),
              by_sctransform        => getValue( $def, "by_sctransform" ),
              resolution            => getValue( $def, "multires_resolution" ),
              reduction             => $reduction,
              species               => getValue( $def, "species" ),
              db_markers_file       => getValue( $def, "markers_file" ),
              curated_markers_file  => getValue( $def, "curated_markers_file", "" ),
              annotate_tcell        => getValue( $def, "annotate_tcell", 0),
              remove_subtype        => getValue( $def, "remove_subtype", ""),
              HLA_panglao5_file     => getValue( $def, "HLA_panglao5_file", "" ),
              tcell_markers_file    => getValue( $def, "tcell_markers_file", ""),
              summary_layer_file => $def->{summary_layer_file},
            },
            parameterSampleFile2 => getValue($def, "subcluster_ignore_gene_files", {}),
            output_file_ext      => ".final.rds",
            output_other_ext  => ".umap.sample_cell.png",
            sh_direct            => 1,
            pbs                  => {
              "nodes"     => "1:ppn=1",
              "walltime"  => "1",
              "mem"       => "40gb"
            },
          };
          push( @$summary, $recluster_task );
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
          parameterFile2_ref => [$seurat_task, ".pca.txt"],
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
            "walltime"  => "1",
            "mem"       => "40gb"
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
              "walltime"  => "1",
              "mem"       => "40gb"
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
          addScMRMA( $config, $def, $summary, $target_dir, $project_name, $cluster_task );
        }

        my $parameterSampleFile5_ref = undef;
        if (getValue( $def, "perform_CHETAH", 0 ) ) {
          my $CHETAH_name= $cluster_task . "_CHETAH";
          addCHETAH( $config, $def, $summary, $target_dir, $project_name, $CHETAH_name, $cluster_task );
          push @report_files, ($CHETAH_name, ".CHETAH.png");
          push @report_names, "chetah_png";
        }

        addGeneTask( $config, $def, $summary, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name );

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

          addGeneTask( $config, $def, $summary, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name );
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
          addAntibody( $config, $def, $summary, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name );
        }

        if(defined $clonotype_4_convert){
          addClonotypeVis($config, $def, $summary, $target_dir, $clonotype_4_convert . "_vis", [ $cluster_task, ".final.rds" ], [ $celltype_task, $celltype_cluster_file ], $clonotype_4_convert);

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
        
        if ( $perform_comparison ) {
          if ( defined $def->{"DE_cluster_pairs"} ) {
            addEdgeRTask( $config, $def, $summary, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name, 1, 0, $DE_by_cell );
          }

          for my $deByOption (@deByOptions) {
            my $DE_by_celltype = $deByOption eq "DE_by_celltype";
            addEdgeRTask( $config, $def, $summary, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name, 0, $DE_by_celltype, $DE_by_cell );
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

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => $individual,
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
