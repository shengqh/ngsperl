#!/usr/bin/perl
package scRNA::Modules;

use strict;
use warnings;
require Exporter;
use CQS::ConfigUtils;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(addEnclone addClonotypeMerge addEncloneToClonotype addArcasHLA addScMRMA addCellRangerCount addCellRangerVdj)] );

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
  push @$tasks, $taskname;
  return($taskname);
}

sub addArcasHLA {
  my ( $config, $def, $tasks, $target_dir, $task_name, $prefix, $source_ref ) = @_;

  $config->{"${prefix}_arcasHLA_1_extract"} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/${prefix}_arcasHLA_1_extract",
    init_command          => "ln -s __FILE__ __NAME__.bam",
    option                => "extract -t 8 -v --log __NAME__.log --paired __NAME__.bam",
    interpretor           => "",
    check_program         => 0,
    program               => "arcasHLA",
    source_ref            => $source_ref,
    source_arg            => "--paired",
    source_join_delimiter => "",
    output_to_same_folder => 0,
    output_arg            => "-o",
    output_to_folder      => 1,
    output_file_prefix    => "",
    output_file_ext       => ".extracted.1.fq.gz",
    output_other_ext      => ".extracted.2.fq.gz",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };

  $config->{"${prefix}_arcasHLA_2_genotype"} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "$target_dir/${prefix}_arcasHLA_2_genotype",
    option                => "genotype -t 8 -v --log __NAME__.log",
    interpretor           => "",
    check_program         => 0,
    program               => "arcasHLA",
    source_ref            => "${prefix}_arcasHLA_1_extract",
    source_arg            => "",
    source_join_delimiter => " ",
    output_to_same_folder => 1,
    output_to_folder      => 1,
    output_arg            => "-o",
    output_file_prefix    => "",
    output_file_ext       => ".genotype.json",
    output_other_ext       => ".genes.json,.alignment.p",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => "40gb"
    },
  };

  $config->{"${prefix}_arcasHLA_3_merge"} = {
    class                 => "CQS::ProgramWrapper",
    perform               => 1,
    target_dir            => "$target_dir/${prefix}_arcasHLA_3_merge",
    option                => "merge -i $target_dir/genotype/result --run $task_name -o .",
    interpretor           => "",
    check_program         => 0,
    program               => "arcasHLA",
    source_ref            => "${prefix}_arcasHLA_2_genotype",
    source_arg            => "-i",
    source_join_delimiter => " ",
    output_to_same_folder => 1,
    output_to_folder      => 1,
    output_arg            => "-o",
    output_file_prefix    => "",
    output_file_ext       => ".genotype.txt",
    sh_direct             => 1,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };

  push (@$tasks, ("${prefix}_arcasHLA_1_extract", "${prefix}_arcasHLA_2_genotype", "${prefix}_arcasHLA_3_merge"));
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
