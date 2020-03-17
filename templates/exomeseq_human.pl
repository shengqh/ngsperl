#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use CQS::PerformExomeSeq;

my $def = {
  task_name => "exomeseq_human_example",

  target_dir => "/scratch/cqs/shengq2/temp/exomeseq_human_example",
  email      => "quanhu.sheng.1\@vumc.org",
  "mail-type"  => "FAIL",
  files      => {
     "MEL194" =>      [ "/data/cqs/macrae_linton_data/2118/2118-JB-1-AGTGTTGC-ATGTAACG_S329_R1_001.fastq.gz", "/data/cqs/macrae_linton_data/2118/2118-JB-1-AGTGTTGC-ATGTAACG_S329_R2_001.fastq.gz" ],
    "AF245" =>      [ "/data/cqs/macrae_linton_data/2118/2118-JB-2-TTACCTGG-GATTCTGA_S330_R1_001.fastq.gz", "/data/cqs/macrae_linton_data/2118/2118-JB-2-TTACCTGG-GATTCTGA_S330_R2_001.fastq.gz" ],
  },
  merge_fastq => 0,

  perform_fastqc => 0,

  is_paired        => 1,
  perform_cutadapt => 0,
  adapter          => "AGATCGGAAGAGC",
  min_read_length  => 30,

  aligner_scatter_count => 5,

  covered_bed                 => "/scratch/cqs_share/references/exomeseq/IDT/Exome-IDT-xGen-hg19-v1-slop50-nochr.bed",
  perform_gatk4_callvariants  => 0,

  gatk4_variant_filter_by_chromosome => 0,

  perform_gatk_callvariants   => 0,
  gatk_callvariants_vqsr_mode => 0,

  filter_variants_by_allele_frequency            => 1,
  filter_variants_by_allele_frequency_percentage => 0.9,
  filter_variants_by_allele_frequency_maf        => 0.3,

  #annotation_genes    => "LDLR APOB PCSK9 LDLRAP1 STAP1 LIPA ABCG5 ABCGB APOE LPA PNPLA5 CH25H INSIG2 SIRT1 LRP2",
  #family_info_file    => "/scratch/cqs/shengq2/macrae_linton/20180913_linton_exomeseq_2118_human_cutadapt/linton_exomeseq_2118.pass.family.txt",
  #family_info_feature => "family",

  perform_cnv_cnMOPs => 0,
  perform_vep        => 0,

  perform_cnv_gatk4_cohort => 0,

  cnv_xhmm_preprocess_intervals => 0,
  perform_cnv_xhmm              => 0,
};

my $config = performExomeSeq_gatk_hg19( $def, 1 );

1;

