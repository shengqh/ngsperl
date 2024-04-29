#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

#define task name, this name will be used as prefix of a few result, such as read count table file name.
my $task_name = "B3436";

#email which will be used for notification if you run through cluster
my $email = "quanhu.sheng\@vanderbilt.edu";

#target dir which will be automatically created and used to save code and result
my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/pipeline_star_deseq2");

#cqs in house software which will be used to generate count table. https://github.com/shengqh/CQS.Tools/
my $cqstools = "/home/shengq2/cqstools/cqstools.exe";

#gene annotation file
my $transcript_gtf = "/scratch/cqs/shengq2/references/gatk/b37/Homo_sapiens.GRCh37.82.MT.gtf";

#gene name to transcription name map file which is generated by "mono cqstools.exe gtf_buildmap".
my $name_map_file = "/scratch/cqs/shengq2/references/gatk/b37/Homo_sapiens.GRCh37.82.MT.map";

#genome sequence file
my $fasta_file = "/scratch/cqs/shengq2/references/gatk/b37/human_g1k_v37.fasta";

#star index
my $star_index  = "/scratch/cqs/shengq2/references/gatk/b37/STAR_index_2.5.2b_ensembl_v82_sjdb99";

#rnaseq jar file
my $rnaseqc_jar = "/scratch/cqs/shengq2/local/bin/RNA-SeQC_v1.1.8.jar";

my $config = {
  general => { task_name => $task_name },
  files   => {
    "JDB-01" => [ "/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-1_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-1_2_sequence.txt.gz" ],
    "JDB-02" => [ "/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-2_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-2_2_sequence.txt.gz" ],
    "JDB-03" => [ "/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-3_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-3_2_sequence.txt.gz" ],
    "JDB-04" => [ "/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-4_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-4_2_sequence.txt.gz" ],
    "JDB-05" => [ "/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-5_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-5_2_sequence.txt.gz" ],
    "JDB-06" => [ "/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-6_1_sequence.txt.gz", "/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-6_2_sequence.txt.gz" ],
  },
  groups => {
    "FED"  => [ "JDB-01", "JDB-02", "JDB-03" ],
    "DMSO" => [ "JDB-04", "JDB-05", "JDB-06" ],
  },
  pairs => {
    "DMSO_vs_FED" => [ "FED", "DMSO" ],
  },
  fastqc => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "files",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqc",
    cqstools   => $cqstools,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  star => {
    class                     => "Alignment::STAR",
    perform                   => 1,
    target_dir                => "${target_dir}/star",
    option                    => "--twopassMode Basic",
    source_ref                => "files",
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
  star_rnaseqc => {
    class          => "QC::RNASeQC",
    perform        => 1,
    target_dir     => "${target_dir}/star_rnaseqc",
    option         => "",
    source_ref     => [ "star", ".bam\$" ],
    jar            => $rnaseqc_jar,
    fasta_file     => $fasta_file,
    transcript_gtf => $transcript_gtf,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },

  star_featurecount => {
    class      => "Count::FeatureCounts",
    perform    => 1,
    target_dir => "${target_dir}/star_featurecount",
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
    target_dir    => "${target_dir}/star_genetable",
    option        => "-k 0 -v 6 -e -o ${task_name}_gene.count",
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
    target_dir           => "${target_dir}/star_genetable_deseq2",
    option               => "",
    source_ref           => "pairs",
    groups_ref           => "groups",
    countfile_ref        => "star_genetable",
    sh_direct            => 1,
    show_DE_gene_cluster => 0,
    pvalue               => 0.05,
    fold_change          => 2.0,
    min_median_read      => 5,
    add_count_one        => 0,
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
    target_dir               => "${target_dir}/star_genetable_correlation",
    rtemplate                => "countTableVisFunctions.R,countTableGroupCorrelation.v2.R",
    output_file              => "parameterSampleFile1",
    output_file_ext          => ".Correlation.png",
    parameterSampleFile1_ref => [ "star_genetable", ".count\$" ],
    parameterSampleFile2_ref => "groups",
    sh_direct                => 1,
    pbs                      => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "10gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => [ "fastqc",         "star",         "star_featurecount" ],
      step2 => [ "fastqc_summary", "star_rnaseqc", "star_genetable", "star_genetable_correlation", "star_genetable_deseq2" ],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);
1;

