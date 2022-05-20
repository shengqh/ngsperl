#!/usr/bin/perl
package scRNA::Modules;

use strict;
use warnings;
require Exporter;
use CQS::ConfigUtils;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(
  get_marker_gene_dict
  addEnclone 
  addClonotypeMerge 
  addEncloneToClonotype 
  addArcasHLA 
  addScMRMA 
  addCHETAH
  addSignac_only
  addSignac
  addCellRangerCount 
  addCellRangerVdj
  addDoubletFinder
  addAntibody
  addMarkerGenes
  addGeneTask
  addEdgeRTask
  addComparison)] );

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

sub addEnclone {
  my ( $config, $def, $tasks, $taskName, $parentDir, $sourceRef ) = @_;

  $config->{$taskName} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$parentDir/$taskName",
    init_command          => '',
    option                => "
dn=`dirname __FILE__`

enclone TCR=\${dn} POUT=__NAME__.csv > __NAME__.log",
    interpretor           => "",
    check_program         => 0,
    program               => "",
    source_ref            => $sourceRef,
    source_arg            => "TCR=",
    source_join_delimiter => " ",
    output_to_same_folder => 1,
    output_to_folder      => 1,
    output_arg            => ">",
    output_file_ext       => ".csv",
    no_docker             => 1,
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "4",
      "mem"       => "10gb"
    },
  };

  push(@$tasks, $taskName);
}

sub addClonotypeMerge {
  my ( $config, $def, $tasks, $target_dir, $taskname, $source_ref ) = @_;

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
    output_no_name           => 1,
    sh_direct                => 1,
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };
  push @$tasks, $taskname;
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
    output_file              => "clonotypes",
    output_file_ext          => ".csv",
    output_no_name           => 1,
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
      "mem"       => "10gb"
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
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push( @$tasks, $task_name );
}

sub addSignac_only {
  my ( $config, $def, $tasks, $target_dir, $project_name, $task_name, $seurat_ref, $reduction ) = @_;

  $config->{$task_name} = {
    class                => "CQS::UniqueR",
    perform              => 1,
    target_dir           => $target_dir . "/" . $task_name,
    rtemplate            => "../scRNA/scRNA_func.r,../scRNA/SignacX_only.r",
    parameterFile1_ref   => $seurat_ref,
    parameterSampleFile1 => {
      species             => getValue( $def, "species" ),
      prefix              => $project_name,
      reduction           => $reduction,
      pca_dims            => getValue( $def, "pca_dims" ),
      bubblemap_file        => $def->{bubblemap_file},
      by_sctransform        => getValue( $def, "by_sctransform" ),
    },
    output_file_ext => ".SignacX.png;.SignacX.rds;.meta.rds",
    sh_direct       => 1,
    pbs             => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
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
      "walltime"  => "1",
      "mem"       => "10gb"
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

  $config->{$task_name} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => "${target_dir}/$task_name",
    program => "cellranger",
    check_program => 0,
    option => " count --disable-ui --id=__NAME__ --transcriptome=$count_reference --fastqs=$fastq_folder --sample=__FILE__ $job_arg $chemistry_arg

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
    output_arg => "",
    output_file_prefix => "/filtered_feature_bc_matrix.h5",
    output_file_ext => "/filtered_feature_bc_matrix.h5",
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
    program => "cellranger",
    check_program => 0,
    option => " vdj --disable-ui --id=__NAME__ --reference=$vdj_reference --fastqs=$fastq_folder --sample=__FILE__ $job_arg $chain_arg

if [[ -s __NAME__/outs ]]; then
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
    },
    output_file_ext      => ".meta.rds",
    output_other_ext  => ".meta.csv,.options.csv",
    sh_direct            => 1,
    pbs                  => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "12",
      "mem"       => "40gb"
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

  my $markergenes_task = $celltype_task . "_" . $marker_name;
  $config->{$markergenes_task} = {
    class              => "CQS::UniqueR",
    perform            => 1,
    target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $markergenes_task,
    rtemplate          => "../scRNA/scRNA_func.r;../scRNA/plot_genes.r",
    parameterFile1_ref => [ $cluster_task, ".final.rds" ],
    parameterFile2     => $marker_file,
    parameterFile3_ref => [ $celltype_task, $celltype_cluster_file ],
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
      "mem"       => "10gb"
    },
  };
  
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
      rtemplate          => "../scRNA/scRNAgenes.r",
      parameterFile1_ref => [ $cluster_task, ".final.rds" ],
      parameterFile3_ref => [ $celltype_task, $celltype_cluster_file ],
      output_file_ext    => ".cluster.csv",
      rCode              => "genes='" . $genes . "'; celltype_name='$celltype_name'; cluster_name='$cluster_name'; dotPlotOnly=$dotPlotOnly;",
      sh_direct          => 1,
      pbs                => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    };
    push( @$summary, $genes_task );
  }
}

sub addEdgeRTask {
  my ( $config, $def, $summary, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name, $bBetweenCluster, $DE_by_celltype, $DE_by_cell ) = @_;
  my $rCodeDic = {
    "pvalue" => getValue( $def, "DE_pvalue" ),
    "useRawPvalue" => getValue( $def, "DE_use_raw_pvalue" ),
    "foldChange" => getValue( $def, "DE_fold_change" ),
    "bBetweenCluster" => $bBetweenCluster,
    "DE_by_cell" => $DE_by_cell,
  };

  my $edgeRtaskname  = $celltype_task . "_edgeR";
  my $groups         = undef;
  my $pairs          = undef;
  my $curClusterName = undef;
  my $curClusterDisplayName = undef;

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
      $rCodeDic->{"filter_minTPM"}=getValue( $def, "DE_by_sample_filter_minTPM" );
      $rCodeDic->{"filter_samplePercentage"}=getValue( $def, "DE_by_sample_filter_cellPercentage" );
      $edgeRtaskname = $edgeRtaskname . "_bySample";
    }

    $groups = getValue( $def, "groups" );
    $pairs  = getValue( $def, "pairs" );
  }
  $rCodeDic->{"cluster_name"} = $curClusterName;

  $config->{$edgeRtaskname} = {
    class                => "CQS::UniqueR",
    perform              => 1,
    target_dir           => $target_dir . "/" . getNextFolderIndex($def) . $edgeRtaskname,
    rtemplate            => "../scRNA/edgeR.r",
    parameterFile1_ref   => [ $cluster_task, ".final.rds" ],
    parameterFile2_ref   => [ $celltype_task, $celltype_cluster_file ],
    parameterSampleFile1 => $groups,
    parameterSampleFile2 => $pairs,
    parameterSampleFile3 => $rCodeDic,
    output_file_ext      => ".edgeR.files.csv",
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
    parameterFile1_ref => [ $cluster_task, ".final.rds" ],
    parameterFile2_ref => [$edgeRtaskname],
    parameterFile3_ref => [ $celltype_task, $celltype_cluster_file ],
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
      parameterFile1_ref => [ $cluster_task, ".final.rds" ],
      parameterFile2_ref => [$edgeRtaskname],
      parameterFile3_ref => [ $celltype_task, $celltype_cluster_file ],
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
    my $gsea_makeReport = getValue($def, "gsea_makeReport", 0);

    my $gseaTaskName = $edgeRtaskname . "_GSEA";

    #my $gseaCategories = "'h.all.v6.1.symbols.gmt','c2.all.v6.1.symbols.gmt','c5.all.v6.1.symbols.gmt','c6.all.v6.1.symbols.gmt','c7.all.v6.1.symbols.gmt'";
    $config->{$gseaTaskName} = {
      class                      => "CQS::UniqueR",
      perform                    => 1,
      target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $gseaTaskName,
      docker_prefix              => "gsea_",
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
  }
  
  return ($edgeRtaskname);
}

sub addComparison {
  my ($config, $def, $summary, $target_dir, $cluster_task, $celltype_task, $celltype_cluster_file, $celltype_name, $cluster_name) = @_;
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

1;
