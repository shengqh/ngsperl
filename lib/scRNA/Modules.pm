#!/usr/bin/perl
package scRNA::Modules;

use strict;
use warnings;
require Exporter;
use CQS::ConfigUtils;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(addEnclone 
  addClonotypeMerge 
  addEncloneToClonotype 
  addArcasHLA 
  addScMRMA 
  addCHETAH
  addSignac
  addCellRangerCount 
  addCellRangerVdj)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub addEnclone {
  my ( $config, $def, $tasks, $taskName, $parentDir, $sourceRef ) = @_;

  $config->{$taskName} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$parentDir/$taskName",
    init_command          => '
dn=`dirname __FILE__`
',
    option                => "TCR=\${dn} POUT=__NAME__.csv > __NAME__.log",
    interpretor           => "",
    check_program         => 0,
    program               => "enclone",
    source_ref            => $sourceRef,
    source_arg            => "TCR=",
    source_join_delimiter => " ",
    output_to_same_folder => 1,
    output_to_folder      => 1,
    output_arg            => ">",
    output_file_ext       => ".csv",
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
    init_command          => "ln -s __FILE__ __NAME__.bam 

if [[ -s __FILE__.bai ]]; then
  ln -s __FILE__.bai __NAME__.bam.bai
fi

",
    option                => "extract -t 8 --log __NAME__.log $ispairend_option -v __NAME__.bam -o .

rm __NAME__.bam  
rm __NAME__.bam.bai
",
    interpretor           => "",
    check_program         => 0,
    program               => "arcasHLA",
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
    output_file_ext       => ".genotype.txt",
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

sub addSignac {
  my ( $config, $def, $tasks, $target_dir, $project_name, $task_name, $seurat_name, $tcell_only, $celltype, $reduction ) = @_;

  $config->{$task_name} = {
    class                => "CQS::UniqueR",
    perform              => 1,
    target_dir           => $target_dir . "/" . $task_name,
    rtemplate            => "../scRNA/scRNA_func.r,../scRNA/SignacX.r",
    parameterFile1_ref   => [ $seurat_name, ".final.rds" ],
    parameterSampleFile1 => {
      species             => getValue( $def, "species" ),
      prefix              => $project_name,
      tcell_only          => $tcell_only,
      reduction           => $reduction,
      pca_dims            => getValue( $def, "pca_dims" ),
    },
    output_file_ext => ".SignacX.png;.SignacX.rds",
    sh_direct       => 1,
    pbs             => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };

  if($tcell_only){
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
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => "40gb"
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

1;
