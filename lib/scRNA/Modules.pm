#!/usr/bin/perl
package scRNA::Modules;

use strict;
use warnings;
require Exporter;
use File::Basename;
use File::Slurp;
use CQS::ConfigUtils;
use Pipeline::PipelineUtils;
use Data::Dumper;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(
  get_gmm_demux_option_map

  get_sct_str

  get_marker_gene_dict
  add_seurat_rawdata
  add_seurat_merge_object
  add_seurat
  add_essential_gene
  add_scRNABatchQC

  add_sctk
  add_remove_doublets

  add_hto_samples_preparation  
  add_hto_gmm_demux
  add_hto
  add_hto_summary

  add_souporcell
  add_souporcell_integration

  add_hto_bam

  add_clonotype_split

  add_individual_qc

  add_individual_qc_tasks

  add_individual_dynamic_qc

  add_gliph2

  add_group_umap
  add_pseudo_count

  add_doublet_check
  add_scDblFinder

  addEnclone 
  addClonotypeMerge 
  addEncloneToClonotype 
  addEncloneToConsensus
  addConsensusToImmunarch
  addArcasHLA 
  addScMRMA 
  addCHETAH

  add_singleR_cell
  add_azimuth

  add_singleR

  add_signacx_only
  addSignac

  add_decontX
  
  add_celltype_validation

  addCellRangerCount 
  addCellRangerVdj
  addCellRangerMulti

  addDoubletFinder_individual
  addDoubletFinder
  addAntibody
  addMarkerGenes
  addGeneTask
  addEdgeRTask
  addComparison
  addDynamicCluster
  addSubDynamicCluster  
  addSubCluster
  addSubClusterChoose
  addClonotypeVis
  addClonotypeDB
  addClonotypeCluster
  add_strelka2

  add_clustree_rmd

  add_bubble_plots
  add_bubble_files

  add_multiome_qc

  add_fragment_cells

  add_cellbender

  add_cellbender_with_expected_cells  
)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub get_marker_gene_dict {
  my ( $def ) = @_;

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

  return($marker_genes);
}

sub add_seurat_rawdata {
  my ($config, $def, $summary, $target_dir, $seurat_rawdata, $hto_ref, $hto_sample_file, $files_def, $doublets_ref, $doublet_column, $qc_task, $decontX_ref) = @_;

  my $cur_decontX_ref = getValue($def, "remove_decontX") ? $decontX_ref : undef;

  $config->{$seurat_rawdata} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $seurat_rawdata,
    rtemplate                => "../scRNA/scRNA_func.r;../scRNA/seurat_rawdata.r",
    rReportTemplate => "../scRNA/seurat_rawdata.rmd;reportFunctions.R",
    rmd_ext => ".rawdata.html",
    run_rmd_independent => 1,
    parameterSampleFile1_ref => $files_def,
    parameterSampleFile2     => {
      #for report
      task_name => getValue($def, "task_name"),
      email => getValue($def, "email"),
      #for analysis
      Mtpattern             => getValue( $def, "Mtpattern" ),
      rRNApattern           => getValue( $def, "rRNApattern" ),
      hemoglobinPattern     => getValue( $def, "hemoglobinPattern" ),
      species               => getValue( $def, "species" ),
      pool_sample           => getValue( $def, "pool_sample" ),
      hto_sample_file       => $hto_sample_file,
      hto_ignore_samples => $def->{HTO_ignore_samples},
      ensembl_gene_map_file => $def->{ensembl_gene_map_file},
      keep_seurat_object => getValue( $def, "keep_seurat_object", 0 ),
      seurat_sample_column => $def->{"seurat_sample_column"},
      doublet_column => $doublet_column,
      gene_map_file => $def->{gene_map_file},
      hto_regex => $def->{hto_regex},
    },
    parameterSampleFile3 => $def->{"pool_sample_groups"},
    parameterSampleFile4_ref => $hto_ref,
    parameterSampleFile5_ref => $doublets_ref,
    parameterSampleFile6_ref => $qc_task,
    parameterSampleFile7_ref => $cur_decontX_ref,
    parameterFile1 => $def->{filter_config_file},
    output_file_ext      => ".rawobj.rds",
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => getValue($def, "seurat_walltime"),
      "mem"       => getValue($def, "seurat_mem"),
    },
  };

  push( @$summary, $seurat_rawdata );
}

sub add_seurat_merge_object {
  my ($config, $def, $summary, $target_dir, $seurat_rawdata_task, $source_ref, $doublets_ref, $doublet_column) = @_;

  #print($source_ref);

  $config->{$seurat_rawdata_task} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $seurat_rawdata_task,
    rtemplate                => "../scRNA/scRNA_func.r;../scRNA/seurat_merge_object.r",
    rReportTemplate          => "../scRNA/seurat_data.rmd;reportFunctions.R",
    run_rmd_independent => 1,
    parameterSampleFile1_ref => $source_ref,
    parameterSampleFile2     => {
      Mtpattern             => getValue( $def, "Mtpattern" ),
      rRNApattern           => getValue( $def, "rRNApattern" ),
      hemoglobinPattern     => getValue( $def, "hemoglobinPattern" ),
      sample_pattern        => getValue( $def, "sample_pattern" ),
      keep_seurat_object => getValue( $def, "keep_seurat_object", 0 ),
    },
    parameterSampleFile3 => $def->{merge_object_config},
    output_file_ext      => ".rawobj.rds",
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => getValue($def, "seurat_walltime"),
      "mem"       => getValue($def, "seurat_mem"),
    },
  };

  push( @$summary, $seurat_rawdata_task );
}

sub add_seurat {
  my ($config, $def, $summary, $target_dir, $seurat_rawdata, $essential_gene_task, $no_doublets, $is_preprocessed, $prefix, $qc_filter_config_file) = @_;

  my $seurat_task;
  my $reduction;

  my @sample_names = keys %{$def->{files}};
  my $nsamples = scalar(@sample_names);
  my $by_integration = $nsamples > 1 ? getValue( $def, "by_integration" ) : 0;
  my $by_sctransform = getValue( $def, "by_sctransform" );
  my $use_sctransform_v2 = getValue( $def, "use_sctransform_v2", 1);
  my $sct_str = $by_sctransform ? ($use_sctransform_v2 ? "_sct2": "_sct"):"";
  my $thread = getValue( $def, "sctransform_thread", 1 );
  my $rmd_ext = $by_sctransform ? ($use_sctransform_v2 ? ".sct2": ".sct"):"";

  my $preprocessing_rscript;
  if($by_integration){
    my $integration_by_harmony = getValue( $def, "integration_by_harmony", 1);
    if($integration_by_harmony){
      $seurat_task = "${prefix}seurat${sct_str}_harmony";
      $preprocessing_rscript = "../scRNA/seurat_harmony.r";
      $reduction = "harmony";
      $rmd_ext = $rmd_ext . ".harmony";
    }else{
      $seurat_task = "${prefix}seurat${sct_str}_integration";
      $preprocessing_rscript = "../scRNA/seurat_integration.r";
      $reduction = "pca";
      $rmd_ext = $rmd_ext . ".integration";
    }
  }else{
    $seurat_task = "${prefix}seurat${sct_str}_merge";
    $preprocessing_rscript = "../scRNA/seurat_merge.r";
    $reduction = "pca";
    $rmd_ext = $rmd_ext . ".merge";
  }
  $rmd_ext = $rmd_ext . ".html";

  $config->{$seurat_task} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $seurat_task,
    rtemplate                => "../scRNA/scRNA_func.r;$preprocessing_rscript",
    rReportTemplate          => "../scRNA/seurat_data.rmd;reportFunctions.R",
    run_rmd_independent => 1,
    rmd_ext => $rmd_ext,
    parameterFile1_ref => [$seurat_rawdata, ".rawobj.rds"],
    parameterFile2_ref => $essential_gene_task,
    parameterFile3 => $qc_filter_config_file,
    parameterSampleFile1     => {
      task_name             => getValue( $def, "task_name" ),
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
      use_sctransform_v2 => $use_sctransform_v2,
      batch_for_integration => getValue( $def, "batch_for_integration" ),
      qc_genes              => getValue( $def, "qc_genes", "" ),
      is_preprocessed => $is_preprocessed,
      thread => $thread,
    },
    parameterSampleFile2 =>  $def->{"batch_for_integration_groups"},
    output_file_ext      => ".final.rds,.qc.1.png,.qc.2.png,.qc.3.png,.qc.4.png,.sample_cell.csv,.final.png",
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=$thread",
      "walltime"  => getValue($def, "seurat_walltime"),
      "mem"       => getValue($def, "seurat_mem"),
    },
  };

  push( @$summary, $seurat_task );
  return($seurat_task, $reduction);
}

sub addEnclone {
  my ( $config, $def, $tasks, $taskName, $parentDir, $sourceRef ) = @_;

  $config->{$taskName} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$parentDir/$taskName",
    init_command          => '',
    option                => "
dn=`dirname __FILE__`

#enclone TCR=\${dn} POUT=__NAME__.pchain2.pcell.csv PCELL PCHAINS=2 PCOLS=barcode,group_id,group_ncells,clonotype_id,clonotype_ncells,exact_subclonotype_id,nchains,cdr3_dna1,cdr3_aa1,const1,v_name1,d_name1,j_name1,u_cell1,r_cell1,const2,cdr3_dna2,cdr3_aa2,v_name2,d_name2,j_name2,u_cell2,r_cell2    > __NAME__.pchain2.pcell.log
enclone TCR=\${dn} POUT=__NAME__.pchain4.pcell.csv PCELL PCHAINS=4 PCOLS=barcode,group_id,group_ncells,clonotype_id,clonotype_ncells,exact_subclonotype_id,nchains,cdr3_dna1,cdr3_aa1,const1,v_name1,d_name1,j_name1,u_cell1,r_cell1,cdr3_dna2,cdr3_aa2,const2,v_name2,d_name2,j_name2,u_cell2,r_cell2,cdr3_dna3,cdr3_aa3,const3,v_name3,d_name3,j_name3,u_cell3,r_cell3,cdr3_dna4,cdr3_aa4,const4,v_name4,d_name4,j_name4,u_cell4,r_cell4    > __NAME__.pchain4.pcell.log
#enclone TCR=\${dn} POUT=__NAME__.pchain2.csv PCHAINS=2 PCOLS=group_id,group_ncells,clonotype_id,clonotype_ncells,exact_subclonotype_id,cdr3_dna1,cdr3_aa1,v_name1,d_name1,j_name1,cdr3_dna2,cdr3_aa2,v_name2,d_name2,j_name2,barcodes    > __NAME__.pchain2.log
enclone TCR=\${dn} POUT=__NAME__.pchain4.csv PCHAINS=4 PCOLS=n,group_id,group_ncells,clonotype_id,clonotype_ncells,exact_subclonotype_id,cdr3_dna1,cdr3_aa1,v_name1,d_name1,j_name1,cdr3_dna2,cdr3_aa2,v_name2,d_name2,j_name2,cdr3_dna3,cdr3_aa3,v_name3,d_name3,j_name3,cdr3_dna4,cdr3_aa4,v_name4,d_name4,j_name4,barcodes    > __NAME__.pchain4.log
",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $sourceRef,
    source_arg            => "TCR=",
    source_join_delimiter => " ",
    output_to_same_folder => 1,
    output_to_folder      => 1,
    output_arg            => ">",
    output_file_ext       => ".pchain4.pcell.csv,.pchain4.csv",
    no_docker             => 1,
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "4",
      "mem"       => "10gb"
    },
  };

  push(@$tasks, $taskName);
  return($taskName);
}

sub addClonotypeMerge {
  my ( $config, $def, $tasks, $target_dir, $taskname, $source_ref ) = @_;

  my $post_command = "";
  if(defined $def->{enclone_vdj_reference_tar_gz}){
    $post_command = "if [[ ! -e vdj_reference ]]; then
  mkdir vdj_reference
fi

tar -xzvf $def->{enclone_vdj_reference_tar_gz} -C vdj_reference
";
  }

  $config->{$taskname} = {
    class                    => "CQS::ProgramWrapper",
    perform                  => 1,
    target_dir               => "${target_dir}/$taskname",
    option                   => "",
    interpretor              => "python3",
    program                  => "../scRNA/clonotype_merge.py",
    source_arg               => "-i",
    source_ref               => $source_ref,
    output_arg               => "-o",
    output_file              => "all_contig_annotations",
    output_file_ext          => ".json",
    output_other_ext         => ".json.cdr3",
    samplename_in_result => 0,
    post_command             => $post_command,
    sh_direct                => 1,
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };
  push @$tasks, $taskname;
  return($taskname);
}

sub addEncloneToClonotype {
  my ( $config, $def, $tasks, $target_dir, $taskname, $source_ref, $cdr3_ref ) = @_;

  $config->{$taskname} = {
    class                    => "CQS::ProgramWrapper",
    perform                  => 1,
    target_dir               => "${target_dir}/$taskname",
    option                   => "",
    interpretor              => "python3",
    program                  => "../scRNA/enclone_to_clonotype.py",
    # source_arg               => "-i",
    # source_ref               => $source_ref,
    parameterFile1_arg => "-i",
    parameterFile1_ref => $source_ref,
    parameterFile2_arg => "-c",
    parameterFile2_ref => $cdr3_ref,
    output_arg               => "-o",
    output_file              => ".clonotype",
    output_file_ext          => ".csv,.sub.csv",
    sh_direct                => 1,
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };
  push(@$tasks, $taskname);
  return($taskname);
}

sub addEncloneToConsensus {
  my ( $config, $def, $tasks, $target_dir, $taskname, $source_ref, $cdr3_ref, $celltype_ref ) = @_;

  $config->{$taskname} = {
    class                    => "CQS::ProgramWrapper",
    perform                  => 1,
    target_dir               => "${target_dir}/$taskname",
    option                   => "",
    interpretor              => "python3",
    program                  => "../scRNA/enclone_to_consensus_annotations_celltype.py",
    parameterFile1_arg => "-i",
    parameterFile1_ref => $source_ref,
    parameterFile2_arg => "-c",
    parameterFile2_ref => $cdr3_ref,
    parameterFile3_arg => "--celltype_file",
    parameterFile3_ref => $celltype_ref,
    output_arg               => "-o",
    output_file              => ".meta.list",
    output_file_ext          => ".meta.list",
    output_to_same_folder    => 1,
    sh_direct                => 1,
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };
  push(@$tasks, $taskname);
  return($taskname);
}

sub addConsensusToImmunarch {
  my ( $config, $def, $tasks, $target_dir, $taskname, $source_ref ) = @_;

  $config->{$taskname} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => "${target_dir}/$taskname",
    docker_prefix            => "immunarch_",
    option                   => "",
    rtemplate                => "../CQS/Functions.Rmd;../scRNA/immunarch.Rmd;../scRNA/immunarch.r",
    use_vanilla              => 0,
    parameterFile1_ref => $source_ref,
    output_file_ext          => ".all.html",
    sh_direct                => 1,
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };
  push(@$tasks, $taskname);
  return($taskname);
}


sub addArcasHLA_extract {
  my ( $config, $def, $tasks, $target_dir, $task_name, $extract_task, $source_ref, $ispairend ) = @_;

  my $ispairend_option = $ispairend ? "" : "--single";
  my $output_file_ext = $ispairend ? ".extracted.1.fq.gz" : ".extracted.fq.gz";
  my $output_other_ext = $ispairend ? ".extracted.2.fq.gz" : undef;

  $config->{$extract_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$extract_task",
    init_command          => "
if [[ ! -s __NAME__.bam ]]; then
  ln -s __FILE__ __NAME__.bam 
fi

if [[ ! -s __NAME__.bam.bai ]]; then
  if [[ -s __FILE__.bai ]]; then
    ln -s __FILE__.bai __NAME__.bam.bai
  else
    echo samtools index __NAME__.bam
    samtools index __NAME__.bam
  fi
fi
",
    option                => "
echo arcasHLA extract -t 8 --log __NAME__.log $ispairend_option -v __NAME__.bam -o .
arcasHLA extract -t 8 --log __NAME__.log $ispairend_option -v __NAME__.bam -o .

rm __NAME__.bam  
rm __NAME__.bam.bai
",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    ignore_samples        => getValue($def, "arcasHLA_ignore_samples", []),
    source_ref            => $source_ref,
    source_arg            => "-v",
    source_join_delimiter => "",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_to_folder      => 1,
    output_file_prefix    => "",
    output_file_ext       => $output_file_ext,
    output_other_ext      => $output_other_ext,
    docker_prefix         => "arcashla_",
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };

  push (@$tasks, $extract_task);
}

sub addArcasHLA_genotype {
  my ( $config, $def, $tasks, $target_dir, $task_name, $genotype_task, $source_ref, $genotype_options ) = @_;

  $config->{$genotype_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/$genotype_task",
    option                => "genotype -t 8 -v --log __NAME__.log",
    interpretor           => "",
    check_program         => 0,
    program               => "arcasHLA",
    source_ref            => $source_ref,
    source_arg            => "",
    source_join_delimiter => " ",
    parameterSampleFile2  => $genotype_options,
    parameterSampleFile2_arg => "",
    parameterSampleFile2_type => "array",
    output_to_same_folder => 1,
    output_to_folder      => 1,
    output_arg            => "-o",
    output_file_prefix    => "",
    output_file_ext       => ".genotype.json",
    output_other_ext       => ".genes.json,.alignment.p",
    docker_prefix         => "arcashla_",
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => "40gb"
    },
  };

  push (@$tasks, $genotype_task);
}

sub addArcasHLA {
  my ( $config, $def, $tasks, $target_dir, $task_name, $prefix, $source_ref, $singleend_ref ) = @_;

  my $ispairend = is_paired_end( $def );

  my $extract_task_1 = "${prefix}arcasHLA_1_extract";
  addArcasHLA_extract($config, $def, $tasks, $target_dir, $task_name, $extract_task_1, $source_ref, $ispairend );

  my $result1 = get_result_file($config, $extract_task_1);
  my $genotype_options = {};
  for my $sample_name (keys %$result1){
    if($ispairend){
      $genotype_options->{$sample_name} = [""];
    }else{
      $genotype_options->{$sample_name} = ["--single"];
    }
  }

  my $extract_task = $extract_task_1;
  if (defined $singleend_ref){
    my $extract_task_2 = "${prefix}arcasHLA_1_extract_singleend";
    addArcasHLA_extract($config, $def, $tasks, $target_dir, $task_name, $extract_task_2, $singleend_ref, 0 );

    my $result2 = get_result_file($config, $extract_task_2);
    for my $sample_name (keys %$result2){
      $genotype_options->{$sample_name} = ["--single"];
    }

    $extract_task = [ $extract_task_1, $extract_task_2 ];
  }

  my $genotype_task = "${prefix}arcasHLA_2_genotype";
  addArcasHLA_genotype($config, $def, $tasks, $target_dir, $task_name, $genotype_task, $extract_task, $genotype_options );

  my $merge_task="${prefix}arcasHLA_3_merge";
  $config->{$merge_task} = {
    class                 => "CQS::ProgramWrapper",
    perform               => 1,
    target_dir            => "$target_dir/$merge_task",
    option                => "merge -i $target_dir/$genotype_task/result --run $task_name -o .",
    interpretor           => "",
    check_program         => 0,
    program               => "arcasHLA",
    source_ref            => $genotype_task,
    source_arg            => "-i",
    source_join_delimiter => " ",
    output_to_same_folder => 1,
    output_to_folder      => 1,
    output_arg            => "-o",
    output_file_prefix    => "",
    output_file_ext       => ".genotypes.tsv",
    docker_prefix         => "arcashla_",
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };

  push (@$tasks, $merge_task);
  return($merge_task);
}

sub addScMRMA {
  my ( $config, $def, $tasks, $target_dir, $task_name, $seurat_name ) = @_;

  my $scMRMA_name = $seurat_name . "_scMRMA";
  $config->{$scMRMA_name} = {
    class                => "CQS::UniqueR",
    perform              => 1,
    target_dir           => $target_dir . "/" . $scMRMA_name,
    rtemplate            => "../scRNA/scMRMA.r",
    parameterFile1_ref   => [ $seurat_name, ".final.rds" ],
    parameterFile2_ref   => [ $seurat_name, ".cluster.normByUpQuantile.csv" ],
    parameterSampleFile2 => {
      species             => getValue( $def, "species" ),
      prefix              => $task_name,
      db                  => getValue( $def, "scMRMA_db", "hybrid" ),
      p                   => getValue( $def, "scMRMA_pvalue", "0.05" ),
    },
    output_file_ext => ".scMRMA.rds",
    sh_direct       => 1,
    pbs             => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => getValue($def, "seurat_mem")
    },
  };
  push( @$tasks, $scMRMA_name );
}

sub addCHETAH {
  my ( $config, $def, $tasks, $target_dir, $project_name, $task_name, $seurat_name ) = @_;

  $config->{$task_name} = {
    class                => "CQS::UniqueR",
    perform              => 1,
    target_dir           => $target_dir . "/" . $task_name,
    rtemplate            => "../scRNA/CHETAH.r",
    parameterFile1_ref   => [ $seurat_name, ".final.rds" ],
    parameterFile2   => getValue($def, "chetah_reference_file"),
    parameterFile3   => getValue($def, "chetah_ribosomal_file"),
    parameterSampleFile1 => {
      species             => getValue( $def, "species" ),
      prefix              => $project_name,
    },
    output_file_ext => ".CHETAH.png;.CHETAH.rds;.CHETAH.csv",
    sh_direct       => 1,
    pbs             => {
      "nodes"     => "1:ppn=1",
      "walltime"  => getValue($def, "seurat_walltime"),
      "mem"       => getValue($def, "seurat_mem"),
    },
  };
  push( @$tasks, $task_name );
}

sub addSignac {
  my ( $config, $def, $tasks, $target_dir, $project_name, $task_name, $seurat_ref, $tcell_only, $celltype, $reduction, $celltype_layer, $is_dynamic ) = @_;

  $config->{$task_name} = {
    class                => "CQS::UniqueR",
    perform              => 1,
    target_dir           => $target_dir . "/" . $task_name,
    rtemplate            => "../scRNA/scRNA_func.r,../scRNA/SignacX.r",
    parameterFile1_ref   => $seurat_ref,
    parameterSampleFile1 => {
      species             => getValue( $def, "species" ),
      prefix              => $project_name,
      tcell_only          => $tcell_only,
      reduction           => $reduction,
      pca_dims            => getValue( $def, "pca_dims" ),
      celltype_layer      => $celltype_layer,
      bubblemap_file        => $def->{bubblemap_file},
      by_sctransform        => getValue( $def, "by_sctransform" ),
    },
    output_file_ext => ".SignacX.png;.SignacX.rds;.meta.rds",
    sh_direct       => 1,
    pbs             => {
      "nodes"     => "1:ppn=1",
      "walltime"  => getValue($def, "seurat_walltime"),
      "mem"       => getValue($def, "seurat_mem"),
    },
  };

  if($is_dynamic){
    $config->{$task_name}{parameterFile2_ref} = [ $celltype, ".meta.rds" ];
  }else{
    $config->{$task_name}{parameterFile2_ref} = [ $celltype, ".cluster.csv" ];
    $config->{$task_name}{parameterFile3_ref} = [ $celltype, ".celltype.csv" ];
  }
  push( @$tasks, $task_name );
}

sub add_signacx_only {
  my ( $config, $def, $tasks, $target_dir, $project_name, $task_name, $obj_ref, $reduction, $by_individual_sample ) = @_;

  if(!defined $by_individual_sample){
    $by_individual_sample = 0;
  }

  my $class =  $by_individual_sample ? "CQS::IndividualR" : "CQS::UniqueR";
  $config->{$task_name} = {
    class                => $class,
    perform              => 1,
    target_dir           => $target_dir . "/" . $task_name,
    rtemplate            => "../scRNA/scRNA_func.r,../scRNA/SignacX_only.r",
    parameterSampleFile1_ref => $obj_ref,
    parameterSampleFile2 => {
      species             => getValue( $def, "species" ),
      prefix              => $project_name,
      reduction           => $reduction,
      pca_dims            => getValue( $def, "pca_dims" ),
      bubblemap_file        => $def->{bubblemap_file},
      by_sctransform        => getValue( $def, "by_sctransform" ),
      SignacX_reference_file => $def->{SignacX_reference_file}
    },
    output_file_ext => ".SignacX.png;.SignacX.rds;.meta.rds",
    sh_direct       => 0,
    pbs             => {
      "nodes"     => "1:ppn=1",
      "walltime"  => getValue($def, "signacx_walltime", "10"),
      "mem"       => getValue($def, "signacx_mem", $by_individual_sample ? "40g" : getValue($def, "seurat_mem")),
    },
  };

  push( @$tasks, $task_name );
}

sub add_singleR {
  my ( $config, $def, $tasks, $target_dir, $singleR_task, $obj_ref, $meta_ref, $cur_options ) = @_;

  $config->{$singleR_task} = {
    class                => "CQS::UniqueR",
    perform              => 1,
    target_dir           => $target_dir . "/" . $singleR_task,
    rtemplate            => "../scRNA/scRNA_func.r,../scRNA/SingleR.r",
    parameterFile1_ref   => $obj_ref,
    parameterFile2_ref   => $meta_ref,
    parameterSampleFile1 => merge_hash_left_precedent($cur_options,  {
      species             => getValue( $def, "species" ),
      pca_dims            => getValue( $def, "pca_dims" ),
      bubblemap_file        => $def->{bubblemap_file},
      by_sctransform        => getValue( $def, "by_sctransform" ),
    }),
    output_file_ext => ".SingleR.png;.SingleR.rds;.meta.rds",
    sh_direct       => 1,
    pbs             => {
      "nodes"     => "1:ppn=1",
      "walltime"  => getValue($def, "SingleR_walltime", "10"),
      "mem"       => getValue($def, "seurat_mem"),
    },
  };

  push( @$tasks, $singleR_task );
}

sub add_singleR_cell {
  my ( $config, $def, $tasks, $target_dir, $singleR_task, $obj_ref, $cur_options, $by_individual_sample ) = @_;

  if(!defined $by_individual_sample){
    $by_individual_sample = 0;
  }

  my $target_folder = $target_dir . "/" . $singleR_task;
  my $class =  $by_individual_sample ? "CQS::IndividualR" : "CQS::UniqueR";
  $config->{$singleR_task} = {
    class                => $class,
    perform              => 1,
    target_dir           => $target_folder,
    init_command => "",
    rtemplate            => "../scRNA/scRNA_func.r,../scRNA/SingleR.r",
    parameterSampleFile1_ref   => $obj_ref,
    parameterSampleFile2 => merge_hash_left_precedent($cur_options,  {
      species             => getValue( $def, "species" ),
      pca_dims            => getValue( $def, "pca_dims" ),
      bubblemap_file        => $def->{bubblemap_file},
      by_sctransform        => getValue( $def, "by_sctransform" ),
    }),
    output_file_ext => ".SingleR.png;.SingleR.rds;.meta.rds",
    post_command => "rm -rf .cache",
    #no_docker => 1,
    sh_direct       => 0,
    pbs             => {
      "nodes"     => "1:ppn=1",
      "walltime"  => getValue($def, "SingleR_walltime", "10"),
      "mem"       => getValue($def, "SingleR_mem", $by_individual_sample ? "40g" : getValue($def, "seurat_mem")),
    },
  };

  push( @$tasks, $singleR_task );
}

sub add_azimuth {
  my ( $config, $def, $tasks, $target_dir, $azimuth_task, $obj_ref, $cur_options, $by_individual_sample ) = @_;

  if(!defined $by_individual_sample){
    $by_individual_sample = 0;
  }

  my $target_folder = $target_dir . "/" . $azimuth_task;
  my $class =  $by_individual_sample ? "CQS::IndividualR" : "CQS::UniqueR";
  $config->{$azimuth_task} = {
    class                => $class,
    perform              => 1,
    target_dir           => $target_folder,
    init_command => "",
    rtemplate            => "../scRNA/scRNA_func.r,../scRNA/azimuth.r",
    parameterSampleFile1_ref   => $obj_ref,
    parameterSampleFile2 => merge_hash_left_precedent($cur_options,  {
      species             => getValue( $def, "species" ),
      pca_dims            => getValue( $def, "pca_dims" ),
      Azimuth_ref => getValue($def, "Azimuth_ref"),
      bubblemap_file        => $def->{bubblemap_file},
      by_sctransform        => getValue( $def, "by_sctransform" ),
    }),
    output_file_ext => ".azimuth.png;.azimuth.rds;.meta.rds",
    post_command => "rm -rf .cache",
    no_docker => 1,
    sh_direct       => 0,
    pbs             => {
      "nodes"     => "1:ppn=1",
      "walltime"  => getValue($def, "Azimuth_walltime", "10"),
      "mem"       => getValue($def, "Azimuth_mem", $by_individual_sample ? "40g" : getValue($def, "seurat_mem")),
    },
  };

  #print(Dumper($config));
  push( @$tasks, $azimuth_task );
}

sub add_decontX {
  my ( $config, $def, $tasks, $target_dir, $decontX_task, $files_ref, $raw_files_ref, $cur_options, $by_individual_sample ) = @_;

  if(!defined $by_individual_sample){
    $by_individual_sample = 0;
  }

  my $target_folder = $target_dir . "/" . $decontX_task;
  my $class =  $by_individual_sample ? "CQS::IndividualR" : "CQS::UniqueR";
  my $init_command = "";

  my $remove_decontX = getValue($def, "remove_decontX");
  my $remove_decontX_by_contamination = getValue($def, "remove_decontX_by_contamination");

  my $output_file_ext = ".decontX.meta.rds;.decontX.counts.rds;.decontX.png";
  if($remove_decontX && $remove_decontX_by_contamination > 0){
    $output_file_ext = $output_file_ext . ";.decontX.after.png;.decontX.filtered.csv";
  }

  $config->{$decontX_task} = {
    class                => $class,
    perform              => 1,
    target_dir           => $target_dir . "/" . $decontX_task,
    init_command => $init_command,
    rtemplate            => "../scRNA/scRNA_func.r,../scRNA/decontX.r",
    parameterSampleFile1_ref   => $files_ref,
    parameterSampleFile2_ref   => $raw_files_ref,
    parameterSampleFile3 => merge_hash_left_precedent($cur_options,  {
      species => getValue( $def, "species" ),
      remove_decontX => $remove_decontX,
      remove_decontX_by_contamination => $remove_decontX_by_contamination,
    }),
    output_file_ext => $output_file_ext,
    #no_docker => 1,
    sh_direct       => 0,
    pbs             => {
      "nodes"     => "1:ppn=1",
      "walltime"  => getValue($def, "decontX_walltime", "10"),
      "mem"       => getValue($def, "decontX_mem", $by_individual_sample ? "40g" : getValue($def, "seurat_mem")),
    },
  };

  push( @$tasks, $decontX_task );

  if($remove_decontX && ($remove_decontX_by_contamination > 0)){
    my $decontX_summary_task = $decontX_task . "_summary";
    $config->{$decontX_summary_task} = {
      class => "CQS::UniqueRmd",
      target_dir           => $target_dir . "/" . $decontX_summary_task,
      report_rmd_file => "../scRNA/decontX_summary.rmd",
      additional_rmd_files => "../scRNA/scRNA_func.r;reportFunctions.R",
      option => "",
      parameterSampleFile1 => {
        outFile => getValue($def, "task_name"),
        remove_decontX_by_contamination => $remove_decontX_by_contamination,
      },
      parameterSampleFile2_ref => $decontX_task,
      suffix => ".decontX",
      output_file_ext => ".decontX.html",
      can_result_be_empty_file => 0,
      sh_direct   => 1,
      pbs => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "2",
        "mem"       => "10gb"
      },
    };

    push( @$tasks, $decontX_summary_task );
  }
}

sub add_celltype_validation {
  my ( 
    $config, 
    $def, 
    $tasks, 
    $target_dir, 
    $task_name, 
    $object_ref, 
    $meta_ref, 
    $call_files_ref, 
    $celltype_column, 
    $rmd_ext,    
    $signacX_ref, 
    $singleR_ref, 
    $sctk_ref,
    $decontX_ref,
    $is_choose ) = @_;
    
  my $doublet_column = getValue($def, "validation_doublet_column", getValue($def, "doublet_column", "doubletFinder_doublet_label_resolution_1.5"));

  #print("is_choose=" . $is_choose . "\n");

  my $rmd_file = $is_choose ? "../scRNA/seurat_choose_validation.rmd" : "../scRNA/seurat_scDynamic_validation.rmd";
  $config->{$task_name} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $task_name,
    rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_scDynamic_validation.r",
    rReportTemplate => "$rmd_file,reportFunctions.R",
    rmd_ext => $rmd_ext,
    run_rmd_independent => 1,
    parameterFile1_ref => $object_ref,
    parameterFile2_ref => $meta_ref,
    parameterFile3_ref => $call_files_ref,
    parameterSampleFile1 => {
      task_name => getValue($def, "task_name"),
      pca_dims => getValue( $def, "pca_dims" ),
      by_sctransform => getValue( $def, "by_sctransform" ),
      doublet_column => $doublet_column,
      celltype_column => $celltype_column,
      rmd_ext => $rmd_ext,
      bubblemap_file => getValue($def, "bubblemap_file"),
      species => getValue( $def, "species" ),
    },
    parameterSampleFile2 => $def->{pool_sample_groups},
    parameterSampleFile3_ref => $sctk_ref,
    parameterSampleFile4_ref => $signacX_ref,
    parameterSampleFile5_ref => $singleR_ref,
    parameterSampleFile6_ref => $decontX_ref,
    output_file_ext      => ".meta.rds",
    output_other_ext  => "",
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "12",
      "mem"       => getValue($def, "seurat_mem") 
    },
  };
  push( @$tasks, $task_name );
}

sub addCellRangerCount {
  my ( $config, $def, $tasks, $target_dir, $task_name, $fastq_folder, $count_source, $count_reference, $jobmode, $chemistry ) = @_;
  
  my $chemistry_arg = "";
  if((defined $chemistry) and ($chemistry ne "")){
    $chemistry_arg = "--chemistry=$chemistry";
  }
  
  my $job_arg = "";
  if((defined $jobmode) and ($jobmode ne "")){
    $job_arg = "--jobmode=$jobmode";
  }

  my $sh_direct = $job_arg =~ /slurm/;

  my $count_fastq_folders = getValue($def, "count_fastq_folders", {});
  my $names = $config->{$count_source};
  for my $name (sort keys %$names){
    if(not defined $count_fastq_folders->{$name}){
      if($fastq_folder eq ""){
        die "Cannot find fastq folder for $name, you need to define count_fastq_folder in the configuration file.";
      }
      $count_fastq_folders->{$name} = $fastq_folder;
    }
  }

  $config->{$task_name} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => "${target_dir}/$task_name",
    docker_prefix => "cellranger_",
    program => "cellranger",
    check_program => 0,
    option => " count --disable-ui --id=__NAME__ --transcriptome=$count_reference --fastqs=__FILE2__ --sample=__FILE__ $job_arg $chemistry_arg

if [[ -s __NAME__/outs ]]; then
  rm -rf __NAME__/SC_RNA_COUNTER_CS
  mkdir __NAME__/log
  mv __NAME__/_* __NAME__/log   
  mv __NAME__/outs/* __NAME__
  rm -rf __NAME__/outs
fi

#__OUTPUT__
",
    source_arg => "",
    source_ref => $count_source,
    parameterSampleFile2 => $count_fastq_folders,
    output_arg => "",
    output_file_prefix => "/filtered_feature_bc_matrix.h5",
    output_file_ext => "/filtered_feature_bc_matrix.h5,/metrics_summary.csv,/web_summary.html",
    output_to_same_folder => 1,
    can_result_be_empty_file => 0,
    sh_direct   => $sh_direct,
    pbs => {
      "nodes"     => "1:ppn=" . getValue($def, "cellranger_count_cpu", 8),
      "walltime"  => getValue($def, "cellranger_count_walltime", 48),
      "mem"       => getValue($def, "cellranger_count_mem", "40gb"),
    },
  };

  push(@$tasks, $task_name);
}

sub addCellRangerVdj {
  my ( $config, $def, $tasks, $target_dir, $task_name, $fastq_folder, $vdj_source, $vdj_reference, $jobmode, $chain ) = @_;
  
  my $chain_arg = "";
  if((defined $chain) and ($chain ne "")){
    $chain_arg = "--chain=$chain";
  }
  
  my $job_arg = "";
  if((defined $jobmode) and ($jobmode ne "")){
    $job_arg = "--jobmode=$jobmode";
  }

  my $sh_direct = $job_arg =~ /slurm/;

  $config->{$task_name} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => "${target_dir}/$task_name",
    docker_prefix => "cellranger_",
    program => "cellranger",
    check_program => 0,
    option => " vdj --disable-ui --id=__NAME__ --reference=$vdj_reference --fastqs=$fastq_folder --sample=__FILE__ $job_arg $chain_arg

status=\$?
if [[ \$status -ne 0 ]]; then
  touch __NAME__.failed
else
  rm -rf __NAME__/SC_VDJ_ASSEMBLER_CS
  mkdir __NAME__/log
  mv __NAME__/_* __NAME__/log   
  mv __NAME__/outs/* __NAME__
  rm -rf __NAME__/outs
fi

#__OUTPUT__
",
    source_arg => "",
    source_ref => $vdj_source,
    output_arg => "",
    output_file_prefix => "/all_contig_annotations.json",
    output_file_ext => "/all_contig_annotations.json",
    output_to_same_folder => 1,
    can_result_be_empty_file => 0,
    sh_direct   => $sh_direct,
    pbs => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => "40gb"
    },
  };

  push(@$tasks, $task_name);
}

sub writeCellRangerMultiConig {
  my ($def, $files_name,  $target_dir, $config_template) = @_;
  my $files = $def->{$files_name};

  my @lines = read_file($config_template, chomp => 1);

  my $result = {};

  for my $fname (sort keys %$files){
    my $file_def = $files->{$fname};
    my $fastqs = getValue($file_def, "fastqs");

    my $csv_file = "${target_dir}/${fname}.config.csv";
    open(my $csv_fh, ">$csv_file") or die "Cannot create $csv_file";
    for my $line (@lines){
      print $csv_fh $line . "\n";
      if ($line =~ /fastq_id/){
        for my $ftype (sort keys %$file_def){
          if ($ftype eq "fastqs"){
            next;
          }
          print $csv_fh  $file_def->{$ftype} . ",${fastqs},${ftype}\n";
        }
        last;
      }
    }
    close($csv_fh);

    $result->{$fname} = $csv_file;
  }

  return($result);
}

sub addCellRangerMulti {
    my ( $config, $def, $tasks, $target_dir, $task_name, $files_name, $config_template, $jobmode ) = @_;
      
    my $job_arg = "";
    if((defined $jobmode) and ($jobmode ne "")){
      $job_arg = "--jobmode=$jobmode";
    }

    my $csv_files = writeCellRangerMultiConig($def, $files_name, $target_dir, $config_template);

    my $sh_direct = $job_arg =~ /slurm/;
    $config->{$task_name} = {
      class => "CQS::ProgramWrapperOneToOne",
      target_dir => "${target_dir}/$task_name",
      docker_prefix => "cellranger_",
      program => "cellranger",
      check_program => 0,
      option => " multi --disable-ui --id=__NAME__ --csv=__FILE__ $job_arg

if [[ -s __NAME__/outs ]]; then
  rm -rf __NAME__/SC_MULTI_CS
  mkdir __NAME__/log
  mv __NAME__/_* __NAME__/log   
  mv __NAME__/outs/* __NAME__
  rm -rf __NAME__/outs
fi

#__OUTPUT__

",
      source_arg => "",
      source => $csv_files,
      output_arg => "",
      output_file_prefix => "__NAME__/per_sample_outs/__NAME__/count/sample_filtered_feature_bc_matrix.h5",
      output_file_ext => "__NAME__/per_sample_outs/__NAME__/count/sample_filtered_feature_bc_matrix.h5,__NAME__/per_sample_outs/__NAME__/metrics_summary.csv,__NAME__/per_sample_outs/__NAME__/web_summary.html,__NAME__/multi/count/raw_feature_bc_matrix.h5",
      output_to_same_folder => 1,
      samplename_in_result => 1,
      can_result_be_empty_file => 0,
      sh_direct   => $sh_direct,
      pbs => {
          "nodes"     => "1:ppn=" . getValue($def, "cellranger_count_cpu", 8),
          "walltime"  => getValue($def, "cellranger_count_walltime", 48),
          "mem"       => getValue($def, "cellranger_count_mem", "40gb"),
      },
    };

    push(@$tasks, $task_name);
}

sub addDoubletFinder_individual {
  my ( $config, $def, $tasks, $target_dir, $task_name, $object_ref ) = @_;
  my $by_sctransform = getValue( $def, "by_sctransform" );
  my $threads = $by_sctransform ? getValue($def, "sctransform_thread", 8) : 1;
  $config->{$task_name} = {
    class                    => "CQS::IndividualR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $task_name,
    rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_doublet_finder.r",
    parameterSampleFile1_ref => $object_ref,
    parameterSampleFile2     => {
      pca_dims              => getValue( $def, "pca_dims" ),
      by_sctransform        => getValue( $def, "by_sctransform" ),
      use_sctransform_v2    => getValue( $def, "use_sctransform_v2" ),
      threads => $threads,
    },
    output_file_ext      => ".meta.rds",
    output_other_ext  => ".meta.csv,.options.csv",
    sh_direct            => 0,
    pbs                  => {
      "nodes"     => "1:ppn=${threads}",
      "walltime"  => "12",
      "mem"       => getValue($def, "seurat_mem") 
    },
  };
  push( @$tasks, $task_name );
}

sub addDoubletFinder {
  my ( $config, $def, $tasks, $target_dir, $task_name, $object_ref, $meta_ref ) = @_;
  $config->{$task_name} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $task_name,
    rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_doublet_finder.r",
    parameterFile1_ref       => $object_ref,
    parameterFile2_ref       => $meta_ref,
    parameterSampleFile1     => {
      pca_dims              => getValue( $def, "pca_dims" ),
      by_sctransform        => getValue( $def, "by_sctransform" ),
      use_sctransform_v2    => getValue( $def, "use_sctransform_v2" ),
    },
    output_file_ext      => ".meta.rds",
    output_other_ext  => ".meta.csv,.options.csv",
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "12",
      "mem"       => getValue($def, "seurat_mem") 
    },
  };
  push( @$tasks, $task_name );
}

sub addAntibody {
  my ( $config, $def, $summary, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name ) = @_;

  my $pattern = getValue( $def, "antibody_pattern" );

  my $taskname = $celltype_task . "_antibody_vis";
  $config->{$taskname} = {
    class              => "CQS::UniqueR",
    perform            => 1,
    target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $taskname,
    rtemplate          => "../scRNA/scRNAantibody.r",
    parameterFile1_ref => [ $cluster_task, ".final.rds" ],
    parameterFile3_ref => [ $celltype_task, $celltype_cluster_file ],
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
  my ( $config, $def, $summary, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name, $marker_name, $marker_file, $samples ) = @_;

  my $markergenes_task = $cluster_task . "_" . $marker_name;
  if(defined $celltype_task){
    $markergenes_task = $celltype_task . "_" . $marker_name;
  }

  $config->{$markergenes_task} = {
    class              => "CQS::UniqueR",
    perform            => 1,
    target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $markergenes_task,
    rtemplate          => "../scRNA/scRNA_func.r;../scRNA/plot_genes.r",
    parameterFile1_ref => [ $cluster_task, ".final.rds" ],
    parameterFile2     => $marker_file,
    output_file_ext    => ".cluster.csv",
    parameterSampleFile1 => {
      celltype_name => $celltype_name,
      cluster_name => $cluster_name,
      by_sctransform => getValue($def, "by_sctransform"),
      samples => $samples
    },
    sh_direct          => 1,
    pbs                => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => getValue($def, "seurat_mem") 
    },
  };

  if(defined $celltype_task){
    $config->{$markergenes_task}{parameterFile3_ref} = [ $celltype_task, $celltype_cluster_file ];
  }
  
  push( @$summary, $markergenes_task );
}

sub addGeneTask {
  my ( $config, $def, $summary, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name ) = @_;

  my $marker_genes = get_marker_gene_dict($def);

  for my $key (sort keys %$marker_genes){
    my $file = $marker_genes->{$key}{file};
    my $samples = $marker_genes->{$key}{samples};

    addMarkerGenes($config, $def, $summary, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name, "genes_" . $key, $file, $samples);
  }

  if ( defined $def->{genes} ) {
    my $dotPlotOnly = getValue($def, "genesDotPlotOnly", "0");
    my $genes = $def->{genes};
    $genes =~ s/\n/;/g;
    $genes =~ s/\s/;/g;
    $genes =~ s/;;/;/g;
    my $genes_task = $cluster_task . "_genes";
    $config->{$genes_task} = {
      class              => "CQS::UniqueR",
      perform            => 1,
      target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $genes_task,
      rtemplate          => "../scRNA/scRNA_func.r;../scRNA/scRNAgenes.r",
      parameterFile1_ref => [ $cluster_task, ".final.rds" ],
      output_file_ext    => ".cluster.csv",
      rCode              => "genes='" . $genes . "'; celltype_name='$celltype_name'; cluster_name='$cluster_name'; dotPlotOnly=$dotPlotOnly;",
      sh_direct          => 1,
      pbs                => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => getValue($def, "seurat_mem") 
      },
    };

    if(defined $celltype_task){
      $config->{$genes_task}{parameterFile3_ref} = [ $celltype_task, $celltype_cluster_file ];
    }
    push( @$summary, $genes_task );
  }
}

sub addEdgeRTask {
  my ( $config, $def, $summary, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name, $bBetweenCluster, $DE_by_celltype, $DE_by_cell, $reduction ) = @_;

  if(!defined $reduction){
    $reduction = "umap";
  }

  my $rCodeDic = {
    "email" => getValue($def, "email"),
    "affiliation" => $def->{"affiliation"},
    "task_name" => getValue( $def, "task_name" ),
    "pvalue" => getValue( $def, "DE_pvalue" ),
    "useRawPvalue" => getValue( $def, "DE_use_raw_pvalue" ),
    "foldChange" => getValue( $def, "DE_fold_change" ),
    "bBetweenCluster" => $bBetweenCluster,
    "DE_by_cell" => $DE_by_cell,
    "covariance_file" => $def->{covariance_file},
    "sample_column" => $def->{sample_column},
    "group_column" => $def->{group_column},
    "reduction" => $reduction,
  };

  my $edgeRtaskname  = defined $celltype_task ? $celltype_task . "_edgeR" : $cluster_task . "_edgeR";
  my $groups         = undef;
  my $pairs          = undef;
  my $curClusterName = undef;
  my $curClusterDisplayName = undef;

  my $edgeRscript = "../scRNA/edgeR.r";
  my $edgeRmd =  "../scRNA/edgeR.rmd";
  my $rmd_ext = ".edgeR_by_cell.html";
  if ($bBetweenCluster) {
    $edgeRtaskname  = $edgeRtaskname . "_betweenCluster_byCell";
    $curClusterName = getValue( $def, "DE_cluster_name" );
    $curClusterDisplayName = getValue( $def, "DE_cluster_display_name", $curClusterName );
    $rCodeDic->{"filter_minTPM"} = getValue( $def, "DE_by_cell_filter_minTPM" );
    $rCodeDic->{"filter_cellPercentage"}=getValue( $def, "DE_by_cell_filter_cellPercentage" );
    $groups = getValue( $def, "DE_cluster_groups", {} );
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
      $rCodeDic->{"filter_minTPM"}=getValue( $def, "DE_by_cell_filter_minTPM" );
      $rCodeDic->{"filter_cellPercentage"}=getValue( $def, "DE_by_cell_filter_cellPercentage" );
      $edgeRtaskname = $edgeRtaskname . "_byCell";
    }
    else {
      $rCodeDic->{"filter_min_cell_per_sample"}=getValue( $def, "DE_by_sample_min_cell_per_sample" );
      $edgeRtaskname = $edgeRtaskname . "_bySample";
      $edgeRscript = "../scRNA/edgeR_pseudo.r";
      $rmd_ext = ".edgeR_by_sample.html";
    }

    $rCodeDic->{DE_cluster_pattern} = getValue( $def, "DE_cluster_pattern", "*" );

    $groups = getValue( $def, "groups" );
    $pairs  = getValue( $def, "pairs" );
  }
  $rCodeDic->{"cluster_name"} = $curClusterName;

  $config->{$edgeRtaskname} = {
    class                => "CQS::UniqueR",
    perform              => 1,
    target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $edgeRtaskname,
    rtemplate            => "../scRNA/scRNA_func.r,${edgeRscript}",
    rReportTemplate      => "$edgeRmd;reportFunctions.R",
    rmd_ext => $rmd_ext,     
    run_rmd_independent => 1,
    parameterSampleFile1 => $groups,
    parameterSampleFile2 => $pairs,
    parameterSampleFile3 => $rCodeDic,
    output_file_ext      => ".edgeR.files.csv",
    no_docker            => getValue($def, "no_docker", 0),
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "24",
      "mem"       => getValue($def, "seurat_mem") 
    },
  };

  if (!defined $config->{$cluster_task} && -e $cluster_task) {
    $config->{$edgeRtaskname}{parameterFile1} = $cluster_task;
  }else{
    $config->{$edgeRtaskname}{parameterFile1_ref} = [ $cluster_task, ".final.rds" ];
    if(defined $celltype_task){
      $config->{$edgeRtaskname}{parameterFile2_ref} = [ $celltype_task, $celltype_cluster_file ];
    }
  }
  push( @$summary, $edgeRtaskname );

  my $vistaskname = $edgeRtaskname . "_vis";
  my $vis_rmd_txt = $rmd_ext;
  $config->{$vistaskname} = {
    class              => "CQS::UniqueR",
    perform            => 1,
    target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $vistaskname,
    rtemplate          => "../scRNA/scRNA_func.r,../scRNA/edgeRvis.r",
    rReportTemplate      => "../scRNA/edgeRvis.rmd;reportFunctions.R",
    rmd_ext => ".edgeR_vis.html",
    run_rmd_independent => 1,
    parameterFile2_ref => [$edgeRtaskname],
    output_file_ext    => ".vis.files.csv",
    parameterSampleFile1 => {
      "email" => getValue($def, "email"),
      "affiliation" => $def->{"affiliation"},
      "task_name" => getValue( $def, "task_name" ),
      "sample_column" => getValue($def, "sample_column", "orig.ident"),
      "cluster_name" => $curClusterName,
      "bBetweenCluster" => $bBetweenCluster,
      "DE_by_cell" => $DE_by_cell,
      "reduction" => $reduction,
    },
    sh_direct          => 1,
    pbs                => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => getValue($def, "seurat_mem") 
    },
  };
  if (-f $cluster_task) {
    $config->{$vistaskname}{parameterFile1_ref} = $cluster_task;
  }else{
    $config->{$vistaskname}{parameterFile1_ref} = [ $cluster_task, ".final.rds" ];
    if(defined $celltype_task){
      $config->{$vistaskname}{parameterFile3_ref} = [ $celltype_task, $celltype_cluster_file ];
    }
  }

  push( @$summary, $vistaskname );

  if($bBetweenCluster) {
    my $vistaskname2 = $edgeRtaskname . "_dotplot";
    $config->{$vistaskname2} = {
      class              => "CQS::UniqueR",
      perform            => 1,
      target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $vistaskname2,
      rtemplate          => "../scRNA/scRNA_func.r;../scRNA/edgeRdotplot.r",
      parameterFile2_ref => [$edgeRtaskname],
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
        "walltime"  => "10",
        "mem"       => getValue($def, "seurat_mem") 
      },
    };
    if (-e $cluster_task) {
      $config->{$vistaskname2}{parameterFile1} = $cluster_task;
    }else{
      $config->{$vistaskname2}{parameterFile1_ref} = [ $cluster_task, ".final.rds" ];
      $config->{$vistaskname2}{parameterFile3_ref} = [ $celltype_task, $celltype_cluster_file ];
    }

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
        "email" => getValue($def, "email"),
        "affiliation" => $def->{"affiliation"},
        "task_name" => getValue( $def, "task_name" ),
        organism         => getValue( $def, "webgestalt_organism" ),
        interestGeneType => $def->{interestGeneType},
        referenceSet     => $def->{referenceSet},
      },
      output_file_ext    => ".WebGestaltR.files.csv",
      rCode              => "",
      sh_direct          => 1,
      pbs                => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "24",
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
      parameterSampleFile1 => {
        "email" => getValue($def, "email"),
        "affiliation" => $def->{"affiliation"},
        "task_name" => getValue( $def, "task_name" ),
      },
      output_file_ext    => ".link.files.csv",
      sh_direct                  => 1,
      rCode                      => "",
      pbs                        => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "10gb"
      },
    };
    push( @$summary, $linkTaskName );
  }

  if ( getValue( $def, "perform_gsea" ) ) {
    my $gsea_jar        = $def->{gsea_jar}        or die "Define gsea_jar at definition first";
    my $gsea_db         = $def->{gsea_db}         or die "Define gsea_db at definition first";
    my $gsea_categories = $def->{gsea_categories} or die "Define gsea_categories at definition first";

    my $gsea_chip = $def->{gsea_chip};
    my $gsea_chip_str = defined $gsea_chip ? "gseaChip='$gsea_chip';" : "";

    my $gsea_makeReport = getValue($def, "gsea_makeReport", 0);

    my $use_mouse_gsea_db = getValue($def, "use_mouse_gsea_db", 0);

    my $gseaTaskName = $edgeRtaskname . "_GSEA_" . ($use_mouse_gsea_db ? "Mm" : "Hs");

    #my $gseaCategories = "'h.all.v6.1.symbols.gmt','c2.all.v6.1.symbols.gmt','c5.all.v6.1.symbols.gmt','c6.all.v6.1.symbols.gmt','c7.all.v6.1.symbols.gmt'";
    if(getValue($def, "perform_gsea_by_comparison", 1)){
      $config->{$gseaTaskName} = {
        class                      => "CQS::IndividualR",
        perform                    => 1,
        target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $gseaTaskName,
        docker_prefix              => "gsea_",
        rtemplate                  => "GSEAPerform.R",
        output_to_result_directory => 1,
        output_file_ext            => ".gsea.files.csv",
        parameterFile1_ref         => [ $edgeRtaskname, ".edgeR.files.csv\$" ],
        parameterSampleFile1         => $def->{pairs},
        parameterSampleFile2         => {
          "email" => getValue($def, "email"),
          "affiliation" => $def->{"affiliation"},
          "task_name" => getValue( $def, "task_name" ),
          rmdformats => "readthedown",
        },
        sh_direct                  => 0,
        rCode                      => "$gsea_chip_str gseaDb='" . $gsea_db . "'; gseaJar='" . $gsea_jar . "'; gseaCategories=c(" . $gsea_categories . "); makeReport=" . $gsea_makeReport . ";",
        pbs                        => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "23",
          "mem"       => "10gb"
        },
      };
      push( @$summary, $gseaTaskName );

      my $gsea_report = $gseaTaskName . "_report";
      $config->{$gsea_report} = {
        class                      => "CQS::UniqueR",
        perform                    => 1,
        target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $gsea_report,
        rtemplate                  => "GSEAReport.R",
        rReportTemplate            => "GSEAReport.Rmd;../Pipeline/Pipeline.R;reportFunctions.R",
        run_rmd_independent => 1,
        rmd_ext => ".gsea.html",
        parameterSampleFile1_ref   => [$gseaTaskName],
        parameterSampleFile1_names => ["gsea"],
        parameterSampleFile2 => {
          "email" => getValue($def, "email"),
          "affiliation" => $def->{"affiliation"},
          "task_name" => getValue( $def, "task_name" ),
        },
        remove_empty_parameter => 1,
        output_ext => "gsea_files.csv",
        samplename_in_result => 0,
        sh_direct                  => 1,
        pbs                        => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary, $gsea_report );      
    }else{
      $config->{$gseaTaskName} = {
        class                      => "CQS::UniqueR",
        perform                    => 1,
        target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $gseaTaskName,
        docker_prefix              => "gsea_",
        rtemplate                  => "GSEAPerform.R",
        output_to_result_directory => 1,
        output_file_ext            => ".gsea.files.csv",
        parameterFile1_ref         => [ $edgeRtaskname, ".edgeR.files.csv\$" ],
        parameterSampleFile2         => {
          "email" => getValue($def, "email"),
          "affiliation" => $def->{"affiliation"},
          "task_name" => getValue( $def, "task_name" ),
          rmdformats => "readthedown",
        },
        sh_direct                  => 1,
        rCode                      => "$gsea_chip_str gseaDb='" . $gsea_db . "'; gseaJar='" . $gsea_jar . "'; gseaCategories=c(" . $gsea_categories . "); makeReport=" . $gsea_makeReport . ";",
        pbs                        => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "23",
          "mem"       => "10gb"
        },
      };
      push( @$summary, $gseaTaskName );

      my $gsea_report = $gseaTaskName . "_report";
      $config->{$gsea_report} = {
        class                      => "CQS::UniqueR",
        perform                    => 1,
        target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $gsea_report,
        rtemplate                  => "GSEAReport.R",
        rReportTemplate            => "GSEAReport.Rmd;../Pipeline/Pipeline.R;reportFunctions.R",
        run_rmd_independent => 1,
        rmd_ext => ".gsea.html",
        parameterSampleFile1_ref   => [$gseaTaskName],
        parameterSampleFile1_names => ["gsea"],
        parameterSampleFile2 => {
          "email" => getValue($def, "email"),
          "affiliation" => $def->{"affiliation"},
          "task_name" => getValue( $def, "task_name" ),
        },
        remove_empty_parameter => 1,
        output_ext => "gsea_files.csv",
        samplename_in_result => 0,
        sh_direct                  => 1,
        pbs                        => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary, $gsea_report );

    }
  }
  
  return ($edgeRtaskname);
}

sub addComparison {
  my ($config, $def, $summary, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name, $reduction) = @_;
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

  #print("perform_comparison=" . $perform_comparison , "\n");
  #print("DE_by_sample=" . $DE_by_sample , "\n");
  
  if ( $perform_comparison ) {
    if ( defined $def->{"DE_cluster_pairs"} ) {
      addEdgeRTask( $config, $def, $summary, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name, 1, 0, $DE_by_cell, $reduction );
    }

    for my $deByOption (@deByOptions) {
      my $DE_by_celltype = $deByOption eq "DE_by_celltype";
      addEdgeRTask( $config, $def, $summary, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name, 0, $DE_by_celltype, $DE_by_cell, $reduction );
    }
  }
}

sub addDynamicCluster {
  my ($config, $def, $summary, $target_dir, $scDynamic_task, $seurat_task, $essential_gene_task, $reduction, $by_individual_sample, $by_column, $by_harmony) = @_;

  my $output_file_ext = $by_individual_sample ? ".celltype_cell_num.csv":".iter_png.csv,.scDynamic.meta.rds";
  my $rmd_ext = $by_individual_sample ? ".dynamic_call_individual.html":".dynamic_call.html";
  my $rReportTemplate = $by_individual_sample ? "../scRNA/seurat_scDynamic_one_layer_one_resolution_summary.rmd;../scRNA/seurat_scDynamic_one_layer_one_resolution.rmd;reportFunctions.R":"../scRNA/seurat_scDynamic_one_layer_one_resolution.rmd;reportFunctions.R";
  my $rscript = $by_harmony ? "../scRNA/seurat_scDynamic_one_layer_one_resolution_harmony.r" : "../scRNA/seurat_scDynamic_one_layer_one_resolution.r";

  $config->{$scDynamic_task} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $scDynamic_task,
    rtemplate                => "../scRNA/scRNA_func.r,$rscript",
    rReportTemplate => $rReportTemplate,
    rmd_ext => $rmd_ext,
    run_rmd_independent => 1,
    parameterFile1_ref => [$seurat_task, ".rds"],
    parameterFile3_ref => $essential_gene_task,
    parameterSampleFile1     => {
      task_name             => getValue( $def, "task_name" ),
      pca_dims              => getValue( $def, "pca_dims" ),
      by_sctransform        => getValue( $def, "by_sctransform" ),
      regress_by_percent_mt => getValue( $def, "regress_by_percent_mt" ),
      species               => getValue( $def, "species" ),
      db_markers_file       => getValue( $def, "markers_file" ),
      curated_markers_file  => getValue( $def, "curated_markers_file", "" ),
      annotate_tcell        => getValue( $def, "annotate_tcell", 0),
      remove_subtype => getValue( $def, "remove_subtype", ""),
      HLA_panglao5_file => getValue( $def, "HLA_panglao5_file", "" ),
      tcell_markers_file => getValue( $def, "tcell_markers_file", ""),
      bubblemap_file => $def->{bubblemap_file},
      bubblemap_width => $def->{bubblemap_width},
      bubblemap_height => $def->{bubblemap_height},
      bubblemap_use_order => getValue($def, "bubblemap_use_order", 0),
      summary_layer_file => $def->{summary_layer_file},
      best_resolution_min_markers => getValue( $def, "best_resolution_min_markers" ),
      dynamic_by_one_resolution => getValue( $def, "dynamic_by_one_resolution", 0.2 ),
      redo_harmony          => getValue( $def, "subcluster_redo_harmony", 0),
      layer                 => getValue( $def, "dynamic_layer", "Layer4"),
      reduction             => $reduction,
      by_individual_sample  => $by_individual_sample,
      by_column => $by_column,
    },
    parameterSampleFile2 => $def->{"subcluster_ignore_gene_files"},
    parameterSampleFile3 => $def->{"dynamic_layer_umap_min_dist"},
    parameterSampleFile4 => getValue($def, "dynamic_combine_cell_types", {}),
    output_file_ext      => $output_file_ext,
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "24",
      "mem"       => getValue($def, "seurat_mem")
    },
  };

  push( @$summary, $scDynamic_task );
}

sub addSubDynamicCluster {
  my ($config, $def, $summary, $target_dir, $scSubDynamic_task, $seurat_task, $meta_ref, $essential_gene_task, $reduction, $by_individual_sample, $by_column) = @_;

  my $output_file_ext = $by_individual_sample ? ".celltype_cell_num.csv":".iter_png.csv,.scDynamic.meta.rds";
  my $rmd_ext = $by_individual_sample ? ".dynamic_call_individual.html":".dynamic_call.html";
  my $rReportTemplate = $by_individual_sample ? "../scRNA/seurat_scDynamic_one_layer_one_resolution_summary.rmd;../scRNA/seurat_scDynamic_one_layer_one_resolution.rmd;reportFunctions.R":"../scRNA/seurat_scDynamic_one_layer_one_resolution.rmd;../scRNA/seurat_scDynamic_one_layer_one_resolution_summary.rmd;reportFunctions.R";

  $config->{$scSubDynamic_task} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $scSubDynamic_task,
    rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_scSubDynamic_one_layer_one_resolution.r",
    rReportTemplate => $rReportTemplate,
    rmd_ext => $rmd_ext,
    run_rmd_independent => 1,
    parameterFile1_ref => [$seurat_task, ".rds"],
    parameterFile2_ref => $meta_ref,
    parameterFile3_ref => $essential_gene_task,
    parameterSampleFile1     => {
      task_name             => getValue( $def, "task_name" ),
      pca_dims              => getValue( $def, "pca_dims" ),
      by_sctransform        => getValue( $def, "by_sctransform" ),
      regress_by_percent_mt => getValue( $def, "regress_by_percent_mt" ),
      species               => getValue( $def, "species" ),
      db_markers_file       => getValue( $def, "markers_file" ),
      curated_markers_file  => getValue( $def, "curated_markers_file", "" ),
      annotate_tcell        => getValue( $def, "annotate_tcell", 0),
      remove_subtype => getValue( $def, "remove_subtype", ""),
      HLA_panglao5_file => getValue( $def, "HLA_panglao5_file", "" ),
      tcell_markers_file => getValue( $def, "tcell_markers_file", ""),
      bubblemap_file => $def->{bubblemap_file},
      bubblemap_width => $def->{bubblemap_width},
      bubblemap_height => $def->{bubblemap_height},
      bubblemap_use_order => getValue($def, "bubblemap_use_order", 0),
      summary_layer_file => $def->{summary_layer_file},
      best_resolution_min_markers => getValue( $def, "best_resolution_min_markers" ),
      dynamic_by_one_resolution => getValue( $def, "sub_dynamic_by_one_resolution", 0.2 ),
      redo_harmony          => getValue( $def, "sub_dynamic_redo_harmony", 0),
      init_layer                 => getValue( $def, "sub_dynamic_init_layer", "layer4"),
      final_layer                 => getValue( $def, "sub_dynamic_final_layer", "layer5"),
      reduction             => $reduction,
      by_individual_sample  => $by_individual_sample,
      by_column => $by_column,
    },
    parameterSampleFile2 => $def->{"subcluster_ignore_gene_files"},
    parameterSampleFile3 => $def->{"dynamic_layer_umap_min_dist"},
    parameterSampleFile4 => getValue($def, "dynamic_combine_cell_types", {}),
    output_file_ext      => $output_file_ext,
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "24",
      "mem"       => getValue($def, "seurat_mem")
    },
  };

  push( @$summary, $scSubDynamic_task );
}

sub addSubCluster {
  my (
    $config, 
    $def, 
    $summary, 
    $target_dir, 
    $subcluster_task, 
    $obj_ref, 
    $meta_ref, 
    $essential_gene_task, 
    $cur_options, 
    $rename_map, 
    $signacX_ref, 
    $singleR_ref,
    $rmd_ext) = @_;

  my $by_integration;
  my $integration_by_harmony;
  my $subcluster_redo_harmony;
  if(defined $def->{"subcluster_by_harmony"}){
    if($def->{"subcluster_by_harmony"}){
      $by_integration = 1;
      $integration_by_harmony = 1;
      if(getValue($def, "subcluster_redo_harmony", 0)){
        $subcluster_redo_harmony = 1;
      }else{
        $subcluster_redo_harmony = !getValue( $def, "integration_by_harmony" );
      }

      if (!($subcluster_task =~ /_rh/)){
        $subcluster_task = $subcluster_task . "_rh";
      }
      $cur_options->{reduction} = "harmony";
    }else{
      $by_integration = getValue( $def, "by_integration" );
      $integration_by_harmony = 0;
      $subcluster_redo_harmony = 0;
    }
  }else{
      $by_integration = getValue( $def, "by_integration" );
      $integration_by_harmony = getValue( $def, "integration_by_harmony" );
      $subcluster_redo_harmony = getValue( $def, "subcluster_redo_harmony", 0 );
  }

  $config->{$subcluster_task} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $subcluster_task,
    rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_celltype_subcluster.v3.r",
    rReportTemplate          => "../scRNA/seurat_celltype_subcluster.v3.rmd;reportFunctions.R",
    run_rmd_independent => 1,
    rmd_ext => $rmd_ext,
    parameterFile1_ref => $obj_ref,
    parameterFile2_ref => $meta_ref,
    parameterFile3_ref => $essential_gene_task,
    parameterSampleFile1    => merge_hash_left_precedent($cur_options, {
      task_name             => getValue( $def, "task_name" ),
      pca_dims              => getValue( $def, "pca_dims" ),
      by_sctransform        => getValue( $def, "by_sctransform" ),
      by_integration        => $by_integration,
      by_harmony            => $integration_by_harmony,
      regress_by_percent_mt => getValue( $def, "regress_by_percent_mt" ),
      species               => getValue( $def, "species" ),
      db_markers_file       => getValue( $def, "markers_file" ),
      curated_markers_file  => getValue( $def, "curated_markers_file", "" ),
      annotate_tcell        => getValue( $def, "annotate_tcell", 0),
      remove_subtype        => "", #use all subtypes
      #remove_subtype        => getValue( $def, "remove_subtype", ""),
      HLA_panglao5_file     => getValue( $def, "HLA_panglao5_file", "" ),
      tcell_markers_file    => getValue( $def, "tcell_markers_file", ""),
      redo_harmony          => $subcluster_redo_harmony,
      bubblemap_file => $def->{bubblemap_file},
      antibody_bubblemap_file => $def->{antibody_bubblemap_file},
      bubblemap_width => $def->{bubblemap_width},
      bubblemap_height => $def->{bubblemap_height},
      bubblemap_use_order   => getValue($def, "bubblemap_use_order", 0),
      summary_layer_file    => $def->{summary_layer_file},
      celltype_layer        => "layer4",
      output_layer          => "cell_type",
      best_resolution_min_markers => getValue( $def, "best_resolution_min_markers" ),
      resolutions => getValue( $def, "subcluster_resolutions", "" ),
    }),
    parameterSampleFile2 => $def->{"subcluster_ignore_gene_files"},
    parameterSampleFile3 => $rename_map,
    parameterSampleFile4_ref => $signacX_ref,
    parameterSampleFile5_ref => $singleR_ref,
    parameterSampleFile6 => $def->{dynamic_bubble_files},
    output_file_ext      => ".meta.rds,.files.csv",
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "24",
      "mem"       => getValue($def, "seurat_mem")
    },
  };

  push( @$summary, $subcluster_task );

  return($subcluster_task);
}

sub addSubClusterChoose {
  my (
    $config, 
    $def, 
    $summary, 
    $target_dir, 
    $choose_task, 
    $obj_ref, 
    $meta_ref, 
    $subcluster_task, 
    $essential_gene_task, 
    $cur_options, 
    $celltype_subclusters_table,
    $rmd_ext) = @_;

  $config->{$choose_task} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $choose_task,
    rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_celltype_subcluster.choose.r",
    
    rReportTemplate          => "../scRNA/seurat_celltype_subcluster.choose.rmd;reportFunctions.R",
    run_rmd_independent => 1,
    rmd_ext => $rmd_ext,

    parameterFile1_ref => $obj_ref,
    parameterFile2_ref => $meta_ref,
    parameterFile3_ref => $essential_gene_task,
    parameterFile4_ref => [$subcluster_task, ".files.csv"],
    parameterSampleFile1     => merge_hash_left_precedent($cur_options, {
      task_name => getValue( $def, "task_name" ),
      pca_dims              => getValue( $def, "pca_dims" ),
      by_sctransform        => getValue( $def, "by_sctransform" ),
      regress_by_percent_mt => getValue( $def, "regress_by_percent_mt" ),
      species               => getValue( $def, "species" ),
      db_markers_file       => getValue( $def, "markers_file" ),
      curated_markers_file  => getValue( $def, "curated_markers_file", "" ),
      annotate_tcell        => getValue( $def, "annotate_tcell", 0),
      HLA_panglao5_file     => getValue( $def, "HLA_panglao5_file", "" ),
      tcell_markers_file    => getValue( $def, "tcell_markers_file", ""),
      bubblemap_file => $def->{bubblemap_file},
      bubblemap_width => $def->{bubblemap_width},
      bubblemap_height => $def->{bubblemap_height},
      bubblemap_use_order   => getValue($def, "bubblemap_use_order", 0),
      summary_layer_file    => $def->{summary_layer_file},
      output_layer          => "cell_type",
    }),
    parameterSampleFile3 => $celltype_subclusters_table,
    output_file_ext      => ".meta.rds",
    output_other_ext  => ".final.rds,.meta.csv,.umap.csv,.umap.png",
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "24",
      "mem"       => getValue($def, "seurat_mem")
    },
  };
  push( @$summary, $choose_task );

}

sub addClonotypeVis {
  my ( $config, $def, $tasks, $target_dir, $taskname, $object_ref, $cell_type_ref, $clonotype_convert ) = @_;

  $config->{$taskname} = {
    class                      => "CQS::UniqueR",
    perform                    => 1,
    target_dir                 => $target_dir . "/" . $taskname,
    rtemplate                  => "../scRNA/scRNA_func.r;../scRNA/clonotype_vis.r",
    output_to_result_directory => 1,
    output_file_ext            => ".clonotype_vis.csv",
    parameterFile1_ref         => $object_ref,
    parameterFile2_ref         => $cell_type_ref,
    parameterSampleFile1_ref         => $clonotype_convert,
    sh_direct                  => 1,
    pbs                        => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => getValue($def, "seurat_mem")
    },
  };
  push(@$tasks, $taskname);
}

sub addClonotypeDB {
  my ( $config, $def, $tasks, $target_dir, $taskname, $clonotype_convert ) = @_;
  $config->{$taskname} = {
    class                      => "CQS::UniqueR",
    perform                    => 1,
    target_dir                 => $target_dir . "/" . $taskname,
    rtemplate                  => "../scRNA/clonotype_db.r",
    output_to_result_directory => 1,
    output_file_ext            => ".clonotype.db.csv,.clonotype.sub.db.csv,",
    parameterSampleFile1_ref         => [ $clonotype_convert ],
    parameterFile2         => getValue($def, "clonotype_McPAS_TCR"),
    parameterFile3         => getValue($def, "clonotype_TBAdb"),
    parameterFile4         => getValue($def, "clonotype_vdjdb"),
    sh_direct                  => 1,
    pbs                        => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "10gb"
    },
  };
  push(@$tasks, $taskname);
  return($taskname);
}

sub addClonotypeCluster {
  my ( $config, $def, $tasks, $target_dir, $taskname, $clonotype_db, $celltype_ref, $output_file_ext ) = @_;
  $config->{$taskname} = {
    class                      => "CQS::UniqueR",
    perform                    => 1,
    target_dir                 => $target_dir . "/" . $taskname,
    rtemplate                  => "../scRNA/clonotype_cluster.r",
    output_to_result_directory => 1,
    output_file_ext            => $output_file_ext,
    parameterSampleFile1_ref   => [ $clonotype_db ],
    parameterFile2_ref         => $celltype_ref,
    sh_direct                  => 1,
    pbs                        => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "10gb"
    },
  };
  push(@$tasks, $taskname);
}

sub add_essential_gene {
  my ($config, $def, $summary, $target_dir) = @_;

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

  return($essential_gene_task);
}

sub add_scRNABatchQC{
  my ($config, $def, $summary, $target_dir) = @_;
  my $task = "scRNABatchQC";
  $config->{$task} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $task,
    rtemplate                => "../scRNA/scRNABatchQC.r",
    parameterSampleFile1_ref => "files",
    rCode                    => "webgestalt_organism='" . getValue( $def, "webgestalt_organism" ) . "'",
    output_file_ext      => ".html",
    sh_direct                => 1,
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => getValue($def, "seurat_mem")
    },
  };
  push( @$summary, $task );
  return($task);
}

sub add_hto_samples_preparation {
  my ($config, $def, $summary, $target_dir, $hto_file_ref, $hto_raw_file_ref) = @_;
  my $preparation_task = "hto_samples_preparation";
  $config->{$preparation_task} = {
    class => "CQS::UniqueR",
    target_dir => "${target_dir}/$preparation_task",
    rtemplate => "../scRNA/split_samples_utils.r,../scRNA/scRNA_func.r,../scRNA/split_samples_preparation.r",
    rReportTemplate => "../scRNA/split_samples_preparation.rmd;reportFunctions.R",
    run_rmd_independent => 1,
    rmd_ext => ".hto_preparation.html",
    option => "",
    parameterSampleFile1_ref => $hto_file_ref,
    parameterSampleFile2 => {
      task_name => getValue($def, "task_name"),
      email => getValue($def, "email"),
      hto_regex => getValue($def, "hto_regex", ""),
      nFeature_cutoff_min => getValue($def, "nFeature_cutoff_min"),
      hto_non_zero_percentage => getValue($def, "hto_non_zero_percentage", 0.2),
      hto_filter_by_exp => getValue($def, "hto_filter_by_exp", 1),
      hto_min_count => getValue($def, "hto_min_count", 2),
    },
    #for cellbender result, we need filtered data to extract HTO data from raw data.
    parameterSampleFile3 => $def->{filtered_files},
    #don't use raw object. It would make the normalized count shifted to right.
    #parameterSampleFile3_ref => $hto_raw_file_ref,
    output_perSample_file => "parameterSampleFile1",
    output_perSample_file_byName => 1,
    output_perSample_file_ext => ".hto.rds;.barcodes.tsv",
    output_to_same_folder => 1,
    can_result_be_empty_file => 0,
    sh_direct   => 1,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "40gb"
    },
  };
  push( @$summary, $preparation_task );
  return($preparation_task);
}

sub get_gmm_demux_option_map {
  my $def = shift;
  my $HTO_samples = getValue($def, "HTO_samples");
  my $result = {};
  for my $sample (sort keys %$HTO_samples){
    my $sample_tags = $HTO_samples->{$sample};
    my $tags = [sort keys %$sample_tags];
    my $tag_str = join(',', @$tags);
    $result->{$sample} = $tag_str;
  }
  return($result)
}

sub add_hto_gmm_demux {
  my ($config, $def, $tasks, $target_dir, $hto_file_ref, $hto_sample_file) = @_;

  my $hto_gmm_task = "hto_samples_gmm_demux";
  my $gmm_demux_option_map = get_gmm_demux_option_map($def);

  $config->{$hto_gmm_task} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => "${target_dir}/$hto_gmm_task",
    interpretor => "",
    program => "",
    check_program => 0,
    option => "
res_dir=\$( dirname __FILE__ )
GMM-demux \$res_dir/__NAME__ __FILE2__ -f __NAME__ -o __NAME__
",
    source_arg => "-f",
    source_ref => $hto_file_ref,
    parameterSampleFile2_arg => "",
    parameterSampleFile2 => $gmm_demux_option_map,
    output_arg => "-o",
    output_file_prefix => "",
    output_file_ext => "__NAME__/GMM_full.csv,__NAME__/GMM_full.config",
    output_to_same_folder => 1,
    can_result_be_empty_file => 0,
    sh_direct => 1,
    no_docker => 0,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "20gb"
    },
  };
  push( @$tasks, $hto_gmm_task );

  my $hto_task = "hto_samples_gmm_demux_merge";

  $config->{$hto_task} = {
    class => "CQS::UniqueR",
    target_dir => "${target_dir}/$hto_gmm_task",
    rtemplate => "../scRNA/split_samples_utils.r,../scRNA/split_samples_gmm_demux_merge.r",
    rReportTemplate => "../scRNA/split_samples_summary.rmd;reportFunctions.R",
    run_rmd_independent => 1,
    rmd_ext => ".gmm_demux.html",
    option => "",
    parameterFile1 => $hto_sample_file,
    parameterSampleFile1_ref => [$hto_gmm_task, "GMM_full.csv"],
    parameterSampleFile2_ref => [$hto_gmm_task, "GMM_full.config"],
    parameterSampleFile3 => {
      task_name => getValue($def, "task_name"),
      email => getValue($def, "email"),
      method => "GMM_Demux",
      hto_sample_file => $hto_sample_file,
      umap_min_dist => getValue($def, "hto_umap_min_dist", 0.3),
      umap_num_neighbors => getValue($def, "hto_umap_num_neighbors", 30),
    },
    parameterSampleFile4 => $def->{"HTO_tags"},
    parameterSampleFile5_ref => $hto_file_ref,
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
  push( @$tasks, $hto_task );
  return($hto_task);
}

sub add_hto {
  my ($config, $def, $summary, $target_dir, $hto_file_ref, $hto_sample_file) = @_;

  my $hto_task;
  my $rmd_ext;
  my $method = "";
  my $r_script = undef;
  my $scDemultiplex_cutoff_startval = undef;
  my $scDemultiplex_init_by = undef;
  my $thread = 1;
  if ( getValue($def, "split_hto_samples_by_scDemultiplex", 0)){
    $r_script = "../scRNA/scRNA_func.r,../scRNA/split_samples_utils.r,../scRNA/split_samples_scDemultiplex.r";
    $rmd_ext = ".scDemultiplex.html";
    $method = "scDemultiplex";
    $scDemultiplex_init_by = getValue($def, "scDemultiplex_init_by", "demuxmix");
    $hto_task = "hto_samples_scDemultiplex_$scDemultiplex_init_by";
    if ($scDemultiplex_init_by eq "cutoff"){
      $scDemultiplex_cutoff_startval = getValue($def, "scDemultiplex_cutoff_startval", 0);
    }
    $thread = 5;
  } elsif ( getValue($def, "split_hto_samples_by_cutoff", 0) ) {
    $r_script = "../scRNA/split_samples_utils.r,../scRNA/split_samples_cutoff_all.r";
    $hto_task = "hto_samples_cutoff";
    $rmd_ext = ".cutoff.html";
    $method = "cutoff";
  } else {
    $r_script = "../scRNA/split_samples_utils.r,../scRNA/split_samples_seurat_all.r";
    $hto_task = "hto_samples_HTODemux";
    $rmd_ext = ".HTODemux.html";
    $method = "HTODemux";
  }

  $config->{$hto_task} = {
    class => "CQS::UniqueR",
    target_dir => "${target_dir}/$hto_task",
    rtemplate => $r_script,
    rReportTemplate => "../scRNA/split_samples_summary.rmd;reportFunctions.R",
    rmd_ext => $rmd_ext,
    run_rmd_independent => 1,
    option => "",
    parameterFile1 => $hto_sample_file,
    parameterSampleFile1_ref => $hto_file_ref,
    parameterSampleFile2 => $def->{split_hto_samples_cutoff_point},
    parameterSampleFile3 => {
      task_name => getValue($def, "task_name"),
      email => getValue($def, "email"),
      method => $method,
      init_by => $scDemultiplex_init_by,
      cutoff_startval => $scDemultiplex_cutoff_startval,
      hto_ignore_exists => getValue($def, "hto_ignore_exists", 0),
      cutoff_file => getValue($def, "cutoff_file", ""),
      scDemultiplex_iteration => getValue($def, "scDemultiplex_iteration", 10),
      umap_min_dist => getValue($def, "hto_umap_min_dist", 0.3),
      umap_num_neighbors => getValue($def, "hto_umap_num_neighbors", 30),
    },
    parameterSampleFile4 => $def->{"HTO_tags"},
    output_perSample_file => "parameterSampleFile1",
    output_perSample_file_byName => 1,
    output_perSample_file_ext => ".HTO.umap.class.png;.HTO.csv;.HTO.data.csv;.HTO.umap.rds",
    sh_direct   => 1,
    pbs => {
      "nodes"     => "1:ppn=$thread",
      "walltime"  => "20",
      "mem"       => "40gb"
    },
  };
  push( @$summary, $hto_task );

  return($hto_task);
}

sub add_hto_summary {
  my ($config, $def, $summary, $target_dir, $hto_ref) = @_;
  my $hto_task = $hto_ref->[0];

  my $hto_summary_task = $hto_task . "_summary";
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
  return($hto_summary_task);
}

sub add_souporcell {
  my ($config, $def, $individual, $target_dir, $preparation_task) = @_;

  my $fasta = getValue($def, "souporcell_fasta_file", getValue($def, "", ""));
  my $fasta_map = getValue($def, "souporcell_fasta_file_map", {});

  my $common_variants = getValue($def, "souporcell_common_variants", "");
  my $common_variants_map = getValue($def, "souporcell_common_variants_map", {});

  my $skip_remap = getValue($def, "souporcell_skip_remap", 0);

  my $hto_souporcell_task = "hto_souporcell" . ($skip_remap ? "_skip_remap" : "_remap");

  my $skip_remap_option = $skip_remap ? "--skip_remap SKIP_REMAP" : "";
  my $souporcell_thread = getValue($def, "souporcell_cpu", "16");

  my $bam_files =  getValue($def, "bam_files");
  my $tag_map = getValue($def, "souporcell_tag_number_map");
  my $souporcell_bam_files = {};
  for my $sample (keys %$tag_map){
    $souporcell_bam_files->{$sample} = $bam_files->{$sample};
    if (!defined $fasta_map->{$sample}){
      if($fasta eq ""){
        die "define souporcell_fasta_file_map with $sample, or define souporcell_fasta_file for all samples";
      }else{
        $fasta_map->{$sample} = $fasta;
      }
    }
    if (!defined $common_variants_map->{$sample}){
      if($common_variants eq ""){
        die "define souporcell_common_variants_map with $sample, or define souporcell_common_variants for all samples";
      }else{
        $common_variants_map->{$sample} = $common_variants;
      }
    }
  }

  $config->{$hto_souporcell_task} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => "${target_dir}/$hto_souporcell_task",
    interpretor => "",
    program => "souporcell_pipeline.py",
    check_program => 0,
    docker_prefix => "souporcell_",
    option => "-i __FILE__ -b __FILE2__ -t $souporcell_thread -o . -k __FILE3__ $skip_remap_option",
    post_command => "
#if [[ -s clusters.tsv ]]; then
#  rm -f souporcell_minimap_tagged_sorted.bam.* clusters_tmp.tsv depth_merged.bed minimap.err retag.err clusters.err
#fi
",
    source_arg => "-i",
    source => $souporcell_bam_files,
    parameterSampleFile2_arg => "-b",
    parameterSampleFile2_ref => [ $preparation_task, ".barcodes.tsv"],
    parameterSampleFile3_arg => "-k",
    parameterSampleFile3 => $tag_map,
    parameterSampleFile4_arg => "-f",
    parameterSampleFile4 => $fasta_map,
    parameterSampleFile5_arg => "--common_variants",
    parameterSampleFile5 => $common_variants_map,
    output_arg => "-o",
    output_file_prefix => "",
    samplename_in_result => 0,
    output_file_ext => "clusters.tsv",
    output_to_same_folder => 0,
    can_result_be_empty_file => 0,
    use_tmp_folder => getValue($def, "souporcell_use_tmp_folder", 0),
    sh_direct   => 0,
    pbs => {
      "nodes"     => "1:ppn=" . $souporcell_thread,
      "walltime"  => getValue($def, "souporcell_walltime", "71"),
      "mem"       => getValue($def, "souporcell_mem", "80gb")
    },
  };
  push( @$individual, $hto_souporcell_task );
  return($hto_souporcell_task);
}

sub add_souporcell_integration {
  my ($config, $def, $summary, $target_dir, $hto_souporcell_task, $hto_ref) = @_;
  my $hto_task = $hto_ref->[0];
  my $hto_integration_task = $hto_task . "_souporcell_integration";
  $config->{$hto_integration_task} = {
    class                     => "CQS::UniqueR",
    target_dir                => "${target_dir}/${hto_integration_task}",
    rtemplate                 => "../scRNA/scRNA_func.r;../scRNA/hto_souporcell_integration.r",
    rReportTemplate           => "../scRNA/hto_souporcell_integration.rmd;reportFunctions.R",
    run_rmd_independent => 1,
    rmd_ext => ".souporcell.html",
    option                    => "",
    parameterSampleFile1_ref  => $hto_souporcell_task,
    parameterSampleFile2_ref  => $hto_ref,
    parameterSampleFile3_ref  => [ $hto_task, ".umap.rds" ],
    parameterSampleFile4      => $def->{souporcell_ignore_cluster},
    parameterSampleFile5      => $def->{HTO_samples},
    output_perSample_file     => "parameterSampleFile1",
    output_perSample_file_byName => 1,
    output_perSample_file_ext => ".HTO.png;.HTO.csv;.meta.rds",
    sh_direct                 => 1,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push( @$summary, $hto_integration_task );
  return($hto_integration_task);
}

sub add_hto_bam {
  my ($config, $def, $individual, $target_dir, $hto_ref) = @_;

  my $bam_files = $def->{bam_files};
  my $HTO_samples = $def->{HTO_samples};

  if ( not defined $bam_files){
    die "Define bam_files for perform_arcasHLA";
  }

  if (not defined $HTO_samples) {
    die "Define HTO_samples for split bam files";
  }

  $config->{HTO_samples} = $HTO_samples;
  $config->{bam_files} = $bam_files;

  my $not_hto_bam = {};
  for my $key (sort keys %$bam_files ){
    if(!defined $HTO_samples->{$key}){
      $not_hto_bam->{$key} = $bam_files->{$key};
    }
  }

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
    docker_prefix => "pysam_",
    sh_direct   => 1,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "24",
      "mem"       => "20gb"
    },
  };
  push( @$individual, $hto_bam_task );

  if(%$not_hto_bam){
    $config->{not_hto_bam} = $not_hto_bam;
    return [$hto_bam_task, "not_hto_bam"];
  }else{
    return($hto_bam_task);
  }
}

sub add_clonotype_split {
  my ($config, $def, $individual, $target_dir, $hto_ref, $clono_key) = @_;
  if (not defined $def->{HTO_samples}) {
    die "Define HTO_samples for split vdj json files";
  }

  $config->{HTO_samples} = $def->{HTO_samples};
  $config->{vdj_json_files} = getValue($def, "vdj_json_files");

  my $split_task = "clonotype". get_next_index($def, $clono_key) . "_split";

  $config->{$split_task} = {
    class => "CQS::ProgramWrapperOneToManyFile",
    target_dir => "${target_dir}/$split_task",
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
  push( @$individual, $split_task );
  return($split_task);
}

sub get_sct_str {
  my $def = shift;
  my $by_sctransform = getValue( $def, "by_sctransform" );
  my $use_sctransform_v2 = getValue( $def, "use_sctransform_v2", 0 );
  my $result = $by_sctransform ? ($use_sctransform_v2 ? "_sct2" : "_sct") : "";
  return($result);
}

sub add_individual_qc {
  my ($config, $def, $summary, $target_dir, $individual_qc_task, $qc_filter_config_file, $qc_files_ref, $perform_split_hto_samples, $hto_ref, $hto_sample_file) = @_;

  if(!defined $qc_filter_config_file){
    $qc_filter_config_file = "";
  }
  
  if(!defined $qc_files_ref){
    if(defined $def->{qc_files}){
      $qc_files_ref = "qc_files";
      $config->{qc_files} = $def->{qc_files};
    }else{
      $qc_files_ref = "files";
    }
  }

  #since singleR and signacX depended on the individual qc v2, I will set v2 as default.
  my $v2 = 1;
  #my $v2 = getValue($def, "individual_qc_v2", 1);
  my $class = $v2 ? "CQS::IndividualR" : "CQS::UniqueR";
  my $output_object = $v2 ? 1 : 0;
  my $rtemplate = "../scRNA/individual_qc.r";

  $config->{$individual_qc_task} = {
    class => $class,
    target_dir => "${target_dir}/${individual_qc_task}",
    rtemplate => $rtemplate,
    rReportTemplate => "../scRNA/individual_qc.Rmd;reportFunctions.R;../scRNA/markerCode_filter.R;../scRNA/scRNA_func.r",
    run_rmd_independent => 1,
    rmd_ext => ".individual_qc.html",
    option => "",
    parameterSampleFile1_ref => $qc_files_ref,
    parameterSampleFile2 => {
      email => $def->{email},
      pca_dims              => getValue( $def, "pca_dims" ),
      by_sctransform        => getValue( $def, "by_sctransform" ),
      use_sctransform_v2 => getValue( $def, "use_sctransform_v2", 0 ),
      regress_by_percent_mt => getValue( $def, "regress_by_percent_mt" ),
      species               => getValue( $def, "species" ),
      db_markers_file       => getValue( $def, "markers_file" ),
      curated_markers_file  => getValue( $def, "curated_markers_file", "" ),
      annotate_tcell        => getValue( $def, "annotate_tcell", 0),
      HLA_panglao5_file     => getValue( $def, "HLA_panglao5_file" ),
      tcell_markers_file    => getValue( $def, "tcell_markers_file", ""),
      bubblemap_file => $def->{bubblemap_file},
      bubblemap_width => $def->{bubblemap_width},
      bubblemap_height => $def->{bubblemap_height},
      Mtpattern             => getValue( $def, "Mtpattern" ),
      rRNApattern           => getValue( $def, "rRNApattern" ),
      Remove_rRNA        => getValue( $def, "Remove_rRNA" ),
      Remove_MtRNA        => getValue( $def, "Remove_MtRNA" ),
      nFeature_cutoff_min   => getValue( $def, "nFeature_cutoff_min" ),
      nFeature_cutoff_max   => getValue( $def, "nFeature_cutoff_max" ),
      nCount_cutoff         => getValue( $def, "nCount_cutoff" ),
      mt_cutoff             => getValue( $def, "mt_cutoff" ),
      resolution            => getValue( $def, "resolution" ),
      pca_dims              => getValue( $def, "pca_dims" ),
      ensembl_gene_map_file => $def->{"ensembl_gene_map_file"},
      output_object => $output_object,
    },
    parameterFile1 => $qc_filter_config_file,
    output_file_ext => "objectlist.rds",
    samplename_in_result => 0,
    can_result_be_empty_file => 0,
    remove_empty_parameter => 1,
    sh_direct => 0,
    pbs => {
      "nodes" => "1:ppn=1",
      "walltime" => "10",
      "mem" => getValue($def, "qc_mem", "80gb"),
    },
  };

  if($perform_split_hto_samples){
    $config->{$individual_qc_task}{parameterSampleFile3_ref} = $hto_ref;
    $config->{$individual_qc_task}{parameterSampleFile2}{hto_sample_file} = $hto_sample_file;
  }
    
  if($qc_filter_config_file ne ""){
    if( ! -e $qc_filter_config_file){
      open(my $qc, '>', $qc_filter_config_file) or die $!;
      print $qc "sample,nFeature_cutoff_min,nFeature_cutoff_max,nCount_cutoff,mt_cutoff,cluster_remove\n";
      my $files = $def->{files};
      for my $fname (sort keys %$files){
        print $qc "$fname," . 
                  getValue( $def, "nFeature_cutoff_min" ) . "," . 
                  getValue( $def, "nFeature_cutoff_max" ) . "," . 
                  getValue( $def, "nCount_cutoff" ) . "," .
                  getValue( $def, "mt_cutoff" ) . ",\n";
      }
      close($qc);
    }
  }
  
  push( @$summary, $individual_qc_task );
}

sub add_individual_dynamic_qc {
  my ($config, $def, $summary, $target_dir, $individual_dynamic_qc_task, $qc_filter_config_file, $qc_files_ref, $essential_gene_task) = @_;

  if(!defined $qc_filter_config_file){
    $qc_filter_config_file = "";
  }
  
  if(!defined $qc_files_ref){
    if(defined $def->{qc_files}){
      $qc_files_ref = "qc_files";
      $config->{qc_files} = $def->{qc_files};
    }else{
      $qc_files_ref = "files";
    }
  }

  my $class = "CQS::IndividualR";
  my $output_object = 1;
  my $rtemplate = "../scRNA/individual_dynamic_qc.r";

  $config->{$individual_dynamic_qc_task} = {
    class => $class,
    target_dir => "${target_dir}/${individual_dynamic_qc_task}",
    rtemplate => $rtemplate,
    rReportTemplate => "../scRNA/individual_dynamic_qc.Rmd;reportFunctions.R;../scRNA/markerCode_filter.R;../scRNA/scRNA_func.r",
    run_rmd_independent => 1,
    rmd_ext => ".individual_dynamic_qc.html",
    option => "",
    parameterSampleFile1_ref => $qc_files_ref,
    parameterSampleFile2 => {
      email => $def->{email},
      task_name             => getValue( $def, "task_name" ),
      pca_dims              => getValue( $def, "pca_dims" ),
      by_sctransform        => getValue( $def, "by_sctransform" ),
      regress_by_percent_mt => getValue( $def, "regress_by_percent_mt" ),
      species               => getValue( $def, "species" ),
      db_markers_file       => getValue( $def, "markers_file" ),
      curated_markers_file  => getValue( $def, "curated_markers_file", "" ),
      annotate_tcell        => getValue( $def, "annotate_tcell", 0),
      remove_subtype => getValue( $def, "remove_subtype", ""),
      HLA_panglao5_file => getValue( $def, "HLA_panglao5_file", "" ),
      tcell_markers_file => getValue( $def, "tcell_markers_file", ""),
      bubblemap_file => $def->{bubblemap_file},
      bubblemap_width => $def->{bubblemap_width},
      bubblemap_height => $def->{bubblemap_height},
      bubblemap_use_order => getValue($def, "bubblemap_use_order", 0),
      summary_layer_file => $def->{summary_layer_file},
      best_resolution_min_markers => getValue( $def, "best_resolution_min_markers" ),
      dynamic_by_one_resolution => getValue( $def, "dynamic_by_one_resolution", 0.2 ),
      redo_harmony          => getValue( $def, "subcluster_redo_harmony", 0),
      layer                 => getValue( $def, "dynamic_layer", "Layer4"),

      Mtpattern             => getValue( $def, "Mtpattern" ),
      rRNApattern           => getValue( $def, "rRNApattern" ),
      Remove_rRNA        => getValue( $def, "Remove_rRNA" ),
      Remove_MtRNA        => getValue( $def, "Remove_MtRNA" ),
      nFeature_cutoff_min   => getValue( $def, "nFeature_cutoff_min" ),
      nFeature_cutoff_max   => getValue( $def, "nFeature_cutoff_max" ),
      nCount_cutoff         => getValue( $def, "nCount_cutoff" ),
      mt_cutoff             => getValue( $def, "mt_cutoff" ),
    },
    parameterSampleFile4 => getValue($def, "dynamic_combine_cell_types", {}),
    parameterFile1 => $qc_filter_config_file,
    parameterFile3_ref => $essential_gene_task,
    output_file_ext => ".obj.rds",
    samplename_in_result => 1,
    can_result_be_empty_file => 0,
    remove_empty_parameter => 1,
    sh_direct => 0,
    pbs => {
      "nodes" => "1:ppn=1",
      "walltime" => "10",
      "mem" => getValue($def, "qc_mem", "80gb"),
    },
  };

  if($qc_filter_config_file ne ""){
    if( ! -e $qc_filter_config_file){
      open(my $qc, '>', $qc_filter_config_file) or die $!;
      print $qc "sample,nFeature_cutoff_min,nFeature_cutoff_max,nCount_cutoff,mt_cutoff,cluster_remove\n";
      my $files = $def->{files};
      for my $fname (sort keys %$files){
        print $qc "$fname," . 
                  getValue( $def, "nFeature_cutoff_min" ) . "," . 
                  getValue( $def, "nFeature_cutoff_max" ) . "," . 
                  getValue( $def, "nCount_cutoff" ) . "," .
                  getValue( $def, "mt_cutoff" ) . ",\n";
      }
      close($qc);
    }
  }
  
  push( @$summary, $individual_dynamic_qc_task );
}

sub add_sctk {
  my ($config, $def, $summary, $target_dir, $sctk_task, $files_ref) = @_;

  $config->{$sctk_task} = {
    class                => "CQS::IndividualR",
    perform              => 1,
    target_dir           => $target_dir . "/" . $sctk_task,
    rtemplate            => "../scRNA/scRNA_func.r,../scRNA/sctk.r",
    parameterSampleFile1_ref   => $files_ref,
    parameterSampleFile2 => {
      nFeature_cutoff_min   => getValue( $def, "nFeature_cutoff_min" ),
      nCount_cutoff         => getValue( $def, "nCount_cutoff" ),
    },
    output_file_ext => ".meta.rds",
    #no_docker => 1,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "80gb"
    },
  };

  push( @$summary, $sctk_task );
}

sub add_sctk_old {
  my ($config, $def, $summary, $target_dir, $sctk_task, $files_ref) = @_;

  $config->{$sctk_task} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => "${target_dir}/$sctk_task",
    option => "
echo -e \"__FILE__\\t__NAME__\"> fileList1.txt

R --vanilla -f sctk.r

#__OUTPUT__
",
    check_program => 0,
    program => "",
    source_ref => $files_ref,
    copy_files => "../scRNA/scRNA_func.r;../scRNA/sctk.r",
    output_to_same_folder => 0,
    output_file_ext => ".meta.rds",
    sh_direct   => 0,
    no_docker => 1,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "80gb"
    },
  };

  push( @$summary, $sctk_task );
}

sub add_remove_doublets {
  my ($config, $def, $summary, $target_dir, $task_name, $files_ref, $doublet_meta_ref, $doublet_column) = @_;

  $config->{$task_name} = {
    class                => "CQS::IndividualR",
    perform              => 1,
    target_dir           => $target_dir . "/" . $task_name,
    rtemplate            => "../scRNA/scRNA_func.r,../scRNA/seurat_remove_doublets.r",
    parameterSampleFile1_ref   => $files_ref,
    parameterSampleFile2 => {
      doublet_column => $doublet_column,
    },
    parameterSampleFile3_ref => $doublet_meta_ref,
    output_file_ext => ".nodoublets.counts.rds",
    #no_docker => 1,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => getValue($def, "remove_doublets_mem", "20gb")
    },
  };

  push( @$summary, $task_name );
}

sub add_gliph2 {
  my ($config, $def, $summary, $target_dir, $meta_ref, $clonotype_convert, $hla_merge) = @_;

  my $prepare_task = (defined $hla_merge) ? "tcr_hla_data" : "tcr_data";
  $config->{$prepare_task} = {
    class                => "CQS::UniqueR",
    perform              => 1,
    target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $prepare_task,
    rtemplate            => "../scRNA/gliph2_prepare.r",
    parameterFile1_ref   => $meta_ref,
    parameterFile2_ref   => $clonotype_convert,
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
      "mem"       => getValue($def, "seurat_mem")
    },
  };
  if( defined $hla_merge ) {
    $config->{$prepare_task}{parameterFile3_ref} = $hla_merge; 
    $config->{$prepare_task}{output_other_ext} = ".hla.txt"; 
  }
  push( @$summary, $prepare_task );

  my $gliph2_config_file = dirname(__FILE__) . "/../scRNA/gliph2_config.txt";

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
      "mem"       => getValue($def, "seurat_mem")
    },
  };
  if( defined $hla_merge ) {
    $config->{$gliph2_task}{parameterFile4_ref} = [$prepare_task, '.hla.txt'], 
  }
  push( @$summary, $gliph2_task );
  return($gliph2_task);
}

sub add_group_umap {
  my ($config, $def, $summary, $target_dir, $group_umap_task, $obj_ref) = @_;

  $config->{$group_umap_task} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $group_umap_task,
    rtemplate                => "../scRNA/scRNA_func.r,../scRNA/seurat_group_umap.r",
    parameterFile1_ref       => $obj_ref,
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

sub add_pseudo_count {
  my ($config, $def, $summary, $target_dir, $task_name, $obj_ref, $group_by) = @_;

  $config->{$task_name} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $task_name,
    rtemplate                => "../scRNA/scRNA_func.r,../scRNA/pseudo_count.r",
    parameterFile1_ref       => $obj_ref,
    parameterSampleFile1     => {
      group_by => $group_by
    },
    output_file_ext      => ".pusedo_count.list.csv",
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "12",
      "mem"       => "40gb"
    },
  };
  push( @$summary, $task_name );
}

sub add_doublet_check {
  my ($config, $def, $summary, $target_dir, $doublet_check_task, $obj_ref, $doublet_finder_task ) = @_;

  $config->{$doublet_check_task} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $doublet_check_task,
    rtemplate                => "countTableVisFunctions.R,../scRNA/scRNA_func.r,../scRNA/seurat_doublet_check.r",
    parameterFile1_ref       => $obj_ref,
    parameterFile3_ref       => [ $doublet_finder_task, ".meta.rds" ],
    parameterFile4_ref       => [ $doublet_finder_task, ".options.csv" ],
    parameterSampleFile1     => {
      cluster_layer           => "seurat_clusters",
      celltype_layer          => "cell_type",
      cluster_celltype_layer  => "seurat_cell_type",
      bubblemap_file => $def->{bubblemap_file},
      bubblemap_width => $def->{bubblemap_width},
      bubblemap_height => $def->{bubblemap_height},
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

sub add_scDblFinder {
  my ($config, $def, $summary, $target_dir, $scDblFinder_task, $h5_ref ) = @_;

  my $script = dirname(__FILE__) . "/scDblFinderStandalone.r";

  $config->{$scDblFinder_task} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => "${target_dir}/$scDblFinder_task",
    program => "",
    check_program => 0,
    option => "
R --vanilla -f $script --args __FILE__ __OUTPUT__

",
    source_arg => "",
    source_ref => $h5_ref,
    output_arg => "",
    output_file_prefix => "",
    output_file_ext => ".scDblFinder.rds",
    output_to_same_folder => 1,
    samplename_in_result => 0,
    output_file_key => 0,
    can_result_be_empty_file => 0,
    sh_direct   => 0,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "12",
      "mem"       => "40gb"
    },
  };
  push( @$summary, $scDblFinder_task );
}

#https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1863-4#availability-of-data-and-materials
#https://github.com/fenglin0/benchmarking_variant_callers/blob/master/callVCF_adjustparameters.sh
sub add_strelka2 {
  my ($config, $def, $summary, $target_dir, $strelka2_task, $bam_ref ) = @_;

  my $referenceFasta = getValue($def, "strelka2_referenceFasta");

  $config->{$strelka2_task} = {
    class                    => "CQS::ProgramWrapperOneToOne",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $strelka2_task,
    option                   => "

if [[ -s __FILE__.bai ]]; then
  configureStrelkaGermlineWorkflow.py --rna --bam __FILE__ --referenceFasta $referenceFasta --runDir .

  ./runWorkflow.py -m local -j __THREAD__
fi

#__OUTPUT__
", 
    source_ref => $bam_ref,
    program             => "",
    check_program        => 0,
    output_arg           => "",
    samplename_in_result => 0,
    output_file_ext      => "results/variants/variants.vcf.gz",
    output_other_ext  => "",
    sh_direct            => 1,
    output_to_same_folder => 0,
    no_docker            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => "40gb"
    },
  };
  push( @$summary, $strelka2_task );

  if(getValue($def, "strelka2_extract_snp", 0)){
    my $snp_chrom = getValue($def, "strelka2_extract_snp_chrom");
    my $snp_position = getValue($def, "strelka2_extract_snp_position");

    my $snp_task = $strelka2_task . "_snp";
    $config->{$snp_task} = {
      class                    => "CQS::ProgramWrapper",
      perform                  => 1,
      target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $snp_task,
      option                   => "--snp_chrom $snp_chrom --snp_position $snp_position",
      source_arg => "-i",
      source_ref => $strelka2_task,
      interpretor          => "python3",
      program             => "../Variants/extractVariant.py",
      check_program        => 1,
      output_arg           => "-o",
      output_file_ext      => ".snp.vcf",
      output_other_ext  => "",
      sh_direct            => 1,
      output_to_same_folder => 0,
      no_docker            => 0,
      pbs                  => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "10gb"
      },
    };
    push( @$summary, $snp_task );
  }

  my $combined_task = $strelka2_task . "_combine";
  $config->{$combined_task} = {
    class                    => "CQS::ProgramWrapper",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $combined_task,
    option                   => "

configureStrelkaGermlineWorkflow.py --rna \\
  --bam __FILE__ \\
  --referenceFasta $referenceFasta --runDir .

./runWorkflow.py -m local -j __THREAD__

#__OUTPUT__
", 
    source_arg => "--bam",
    source_type => "array",
    source_join_delimiter => " \\\n  --bam ",
    source_ref => $bam_ref,
    program             => "",
    check_program        => 0,
    output_arg           => "",
    output_file_ext      => ".doublet_perc.png",
    output_other_ext  => "",
    sh_direct            => 1,
    output_to_same_folder => 0,
    no_docker            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "48",
      "mem"       => "60gb"
    },
  };
  push( @$summary, $combined_task );
}

sub add_clustree_rmd {
  my ($config, $def, $summary, $target_dir, $clustree_task, $individual_scDynamic_task, $scDynamic_task) = @_;
  $config->{$clustree_task} = {
    class => "CQS::UniqueRmd",
    target_dir => "${target_dir}/$clustree_task",
    report_rmd_file => "../scRNA/clustree.rmd",
    additional_rmd_files => "../scRNA/scRNA_func.r;reportFunctions.R",
    option => "",
    parameterSampleFile1 => {
      outFile => getValue($def, "task_name")
    },
    parameterSampleFile2_ref => [$scDynamic_task, ".scDynamic.meta.rds"],
    parameterSampleFile3_ref => [$individual_scDynamic_task, ".celltype_cell_num.csv"],
    output_file_ext => "_ur.html",
    can_result_be_empty_file => 0,
    sh_direct   => 1,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "20gb"
    },
  };
  push( @$summary, $clustree_task );
}

sub add_bubble_files {
  my ($config, $def, $summary, $target_dir, $bubble_task, $choose_task, $meta_ref, $celltype_name, $cluster_name, $rmd_ext) = @_;
  my $p2key = defined($meta_ref) ? ((-e $meta_ref) ? "parameterFile2" : "parameterFile2_ref") : "parameterFile2";
  my $suffix = $rmd_ext;
  $suffix =~ s/.html//g;

  $config->{$bubble_task} = {
    class                => "CQS::UniqueR",
    perform              => 1,
    target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $bubble_task,
    rtemplate            => "../scRNA/scRNA_func.r,../scRNA/seurat_bubblemap_multi_slim.r",
    rReportTemplate           => "../scRNA/seurat_bubblemap_multi_slim.rmd;reportFunctions.R",
    rmd_ext => $rmd_ext,
    run_rmd_independent => 1,
    parameterFile1_ref   => [ $choose_task, ".final.rds" ],
    $p2key  => $meta_ref,
    parameterSampleFile1 => $def->{bubble_files},
    parameterSampleFile2 => {
      task_name => getValue($def, "task_name"),
      cluster_name => $cluster_name,
      celltype_name => $celltype_name,
      suffix => $suffix
    },
    output_file_ext      => $rmd_ext,
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "6",
      "mem"       => "40gb"
    },
  };

  push( @$summary, $bubble_task );
}

sub add_bubble_plots {
  my ($config, $def, $summary, $target_dir, $bubble_task, $choose_task, $meta_ref, $celltype_name, $cluster_name, $rmd_ext) = @_;
  my $p2key = defined($meta_ref) ? ((-e $meta_ref) ? "parameterFile2" : "parameterFile2_ref") : "parameterFile2";
  $config->{$bubble_task} = {
    class                => "CQS::UniqueR",
    perform              => 1,
    target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $bubble_task,
    rtemplate            => "../scRNA/scRNA_func.r,../scRNA/seurat_bubblemap_multi.r",
    rReportTemplate           => "../scRNA/seurat_bubblemap_multi.rmd;reportFunctions.R",
    rmd_ext => $rmd_ext,
    run_rmd_independent => 1,
    parameterFile1_ref   => [ $choose_task, ".final.rds" ],
    $p2key  => $meta_ref,
    parameterSampleFile1 => $def->{bubble_plots},
    parameterSampleFile2 => {
      task_name => getValue($def, "task_name"),
      cluster_name => $cluster_name,
      celltype_name => $celltype_name 
    },
    output_file_ext      => ".cell_type.txt",
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "6",
      "mem"       => "40gb"
    },
  };

  push( @$summary, $bubble_task );
}

sub add_individual_qc_tasks{
  my ($config, $def, $summary, $target_dir, $project_name, $prefix, $filter_config_file, $files_def, $decontX_ref, $sctk_ref, $doublet_finder_ref) = @_;

  my $sct_str = get_sct_str($def);
  my $raw_individual_qc_task = "${prefix}raw_qc${sct_str}";

  add_individual_qc($config, $def, $summary, $target_dir, $raw_individual_qc_task, $filter_config_file, $files_def, 0, undef, undef);

  my $reduction = "pca";

  my $signacX_ref = undef;
  if (getValue( $def, "perform_SignacX", 0 ) ) {
    my $signacX_task = $raw_individual_qc_task . "_SignacX";
    add_signacx_only( $config, $def, $summary, $target_dir, $project_name, $signacX_task, $raw_individual_qc_task, $reduction, 1);
    $signacX_ref = [ $signacX_task, ".meta.rds" ];
  }

  my $singleR_ref = undef;
  if (getValue( $def, "perform_SingleR", 0 ) ) {
    my $singleR_task = $raw_individual_qc_task . "_SingleR";
    my $cur_options = {
      task_name => $def->{task_name},
      reduction => $reduction, 
    };
    add_singleR_cell( $config, $def, $summary, $target_dir, $singleR_task, $raw_individual_qc_task, $cur_options, 1 );
    $singleR_ref = [ $singleR_task, ".meta.rds" ];
  }

  my $azimuth_ref = undef;
  if (getValue( $def, "perform_Azimuth", 0 ) ) {
    my $azimuth_task = $raw_individual_qc_task . "_Azimuth";
    my $cur_options = {
      task_name => $def->{task_name},
      reduction => $reduction, 
    };
    add_azimuth( $config, $def, $summary, $target_dir, $azimuth_task, $raw_individual_qc_task, $cur_options, 1);
    $azimuth_ref = [ $azimuth_task, ".meta.rds" ];
    #print($config->{$azimuth_task});
  }

  my $qc_report_task = $raw_individual_qc_task . "_report";
  $config->{$qc_report_task} = {
    class => "CQS::UniqueR",
    perform => 1,
    target_dir => $target_dir . "/" . $qc_report_task,
    rtemplate => "../scRNA/scRNA_func.r,../scRNA/individual_qc_report.r",
    rReportTemplate => "../scRNA/individual_qc_report.Rmd;reportFunctions.R",
    run_rmd_independent => 1,
    rmd_ext => ".${prefix}qc.html",
    parameterSampleFile1 => {
      species => getValue( $def, "species" ),
      prefix => $project_name,
      reduction => $reduction,
      bubblemap_file => $def->{bubblemap_file},
      by_sctransform => getValue( $def, "by_sctransform", 0 ),
      use_sctransform_v2 => getValue( $def, "use_sctransform_v2", 0 ),
      doublet_column => getValue($def, "validation_doublet_column", getValue($def, "doublet_column", "doubletFinder_doublet_label_resolution_1.5")),
    },
    parameterSampleFile2_ref => $raw_individual_qc_task,
    parameterSampleFile3_ref => $sctk_ref,
    parameterSampleFile4_ref => $signacX_ref,
    parameterSampleFile5_ref => $singleR_ref,
    parameterSampleFile6_ref => $decontX_ref,
    parameterSampleFile7_ref => $doublet_finder_ref,
    parameterSampleFile8_ref => $azimuth_ref,
    output_file_ext => ".${prefix}qc.html",
    sh_direct       => 1,
    pbs             => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => getValue($def, "seurat_mem", "40gb") 
    },
  };

  push(@$summary, $qc_report_task);

  return($raw_individual_qc_task, $signacX_ref, $singleR_ref, $qc_report_task);
}

sub add_multiome_qc {
  my ($config, $def, $summary, $target_dir, $multiome_qc_task, $qc_filter_config_file, $fragment_cells_task) = @_;

  if(!defined $qc_filter_config_file){
    $qc_filter_config_file = "";
  }
  
  my $qc_files_ref = "files";
  $config->{atac_files} = $def->{atac_files};
  $config->{raw_files}= $def->{raw_files};

  $config->{$multiome_qc_task} = {
    class => "CQS::IndividualR",
    target_dir => "${target_dir}/${multiome_qc_task}",
    rtemplate => "../scRNA/multiome_qc.r",
    rReportTemplate => "../scRNA/multiome_qc.Rmd;reportFunctions.R;../scRNA/scRNA_func.r",
    run_rmd_independent => 1,
    rmd_ext => ".multiome_qc.html",
    option => "",
    parameterSampleFile1_ref => $qc_files_ref,
    parameterSampleFile2 => {
      email => $def->{email},
      pca_dims              => getValue( $def, "pca_dims" ),
      by_sctransform        => getValue( $def, "by_sctransform" ),
      use_sctransform_v2 => getValue( $def, "use_sctransform_v2", 0 ),
      regress_by_percent_mt => getValue( $def, "regress_by_percent_mt" ),
      species               => getValue( $def, "species" ),
      db_markers_file       => getValue( $def, "markers_file" ),
      curated_markers_file  => getValue( $def, "curated_markers_file", "" ),
      annotate_tcell        => getValue( $def, "annotate_tcell", 0),
      HLA_panglao5_file     => getValue( $def, "HLA_panglao5_file" ),
      tcell_markers_file    => getValue( $def, "tcell_markers_file", ""),
      bubblemap_file => $def->{bubblemap_file},
      bubblemap_width => $def->{bubblemap_width},
      bubblemap_height => $def->{bubblemap_height},
      Mtpattern             => getValue( $def, "Mtpattern" ),
      rRNApattern           => getValue( $def, "rRNApattern" ),
      Remove_rRNA        => getValue( $def, "Remove_rRNA" ),
      Remove_MtRNA        => getValue( $def, "Remove_MtRNA" ),
      nFeature_cutoff_min   => getValue( $def, "nFeature_cutoff_min" ),
      nFeature_cutoff_max   => getValue( $def, "nFeature_cutoff_max" ),
      nCount_cutoff         => getValue( $def, "nCount_cutoff" ),
      mt_cutoff             => getValue( $def, "mt_cutoff" ),
      resolution            => getValue( $def, "resolution" ),
      pca_dims              => getValue( $def, "pca_dims" ),
      ensembl_gene_map_file => $def->{"ensembl_gene_map_file"},

      nCount_ATAC_min => getValue($def, "nCount_ATAC_min", 1000),
      nCount_ATAC_max => getValue($def, "nCount_ATAC_max", 100000),
      max_nucleosome_signal => getValue($def, "max_nucleosome_signal", 2),
      min_TSS_enrichment => getValue($def, "min_TSS_enrichment", 1),
      macs_path => $def->{"macs_path"},
    },
    parameterSampleFile3_ref => "atac_files",
    parameterSampleFile4_ref => "raw_files",
    parameterSampleFile5_ref => $fragment_cells_task,
    parameterFile1 => "",
    parameterFile2 => getValue($def, "signac_annotation_file", ""),
    output_file_ext => ".obj.rds",
    samplename_in_result => 1,
    can_result_be_empty_file => 0,
    remove_empty_parameter => 1,
    sh_direct => 0,
    no_docker => 1,
    pbs => {
      "nodes" => "1:ppn=1",
      "walltime" => "10",
      "mem" => getValue($def, "qc_mem", "80gb"),
    },
  };
}

sub add_fragment_cells {
  my ($config, $def, $tasks, $target_dir, $fragment_cells_task, $fragment_ref ) = @_;

  $config->{$fragment_cells_task} = {
    class                    => "CQS::ProgramWrapperOneToOne",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $fragment_cells_task,
    option                   => "",
    source_arg => "-i", 
    source_ref => $fragment_ref,
    interpretor => "python3",
    program             => "../scRNA/fragment_cells.py",
    check_program        => 1,
    output_arg           => "-o",
    samplename_in_result => 1,
    output_file_prefix => ".cells.txt",
    output_file_ext      => ".cells.txt",
    sh_direct            => 0,
    output_to_same_folder => 1,
    no_docker            => 0,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    },
  };
  push( @$tasks, $fragment_cells_task );
}

sub add_cellbender {
  my ($config, $def, $tasks, $target_dir, $cellbender_task, $files_def ) = @_;

  my $cellbender_cpu = getValue($def, "cellbender_cpu", 12);
  $config->{$cellbender_task} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => "${target_dir}/$cellbender_task",
    program => "",
    check_program => 0,
    option => "
cellbender remove-background --input __FILE__ --output __NAME__.cellbender.h5 --expected-cells __FILE2__ --checkpoint-mins 100000 --cpu-threads $cellbender_cpu

rm -f ckpt.tar.gz

#__OUTPUT__
",
    docker_prefix => "cellbender_",
    parameterSampleFile1_ref => $files_def,
    parameterSampleFile2 => getValue($def, "cellbender_expect_count"),
    output_to_same_folder => 0,
    output_file_ext => ".cellbender_filtered.h5",
    pbs => {
      "nodes"    => "1:ppn=$cellbender_cpu",
      "walltime" => "48",
      "mem"      => "40gb"
    }
  };
  push(@$tasks, $cellbender_task);
}


sub add_cellbender_with_expected_cells {
  my ($config, $def, $tasks, $target_dir, $cellbender_task, $files_def ) = @_;

  my $cellbender_cpu = getValue($def, "cellbender_cpu", 12);
  $config->{$cellbender_task} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => "${target_dir}/$cellbender_task",
    program => "",
    check_program => 0,
    option => "
LOWCOUNTHR=5
NEXP=`cut -d ',' -f4  __FILE2__ | tail -n1`
TOTAL=''
if ((\$NEXP > 20000)) 
then
  NTOT=\$((NEXP+10000))
  echo \"Modifying presets: expecting more than 20k cells (\$NEXP), total number of droplets is \$NTOT..\"
  TOTAL=\"--total-droplets-included \$NTOT\"
else 
  echo \"Standard presets: expected number of cells is \$NEXP..\"
fi

cellbender remove-background --input __FILE__ --output __NAME__.cellbender.h5 --expected-cells \$NEXP \$TOTAL --checkpoint-mins 100000 --cpu-threads $cellbender_cpu

rm -f ckpt.tar.gz

#__OUTPUT__
",
    docker_prefix => "cellbender_",
    parameterSampleFile1_ref => $files_def,
    parameterSampleFile2 => getValue($def, "summary_files"),
    output_to_same_folder => 1,
    output_file_ext => ".cellbender_filtered.h5",
    pbs => {
      "nodes"    => "1:ppn=$cellbender_cpu",
      "walltime" => "48",
      "mem"      => "40gb"
    }
  };
  push(@$tasks, $cellbender_task);
}

1;

