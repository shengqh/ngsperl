#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::ExomeSeq;

my $def = {
  task_name  => "Adenoma",
  target_dir => create_directory_or_die("/scratch/cqs/shengq1/test/20161013_liuqi_gene_panel"),
  email      => "quanhu.sheng.1\@vanderbilt.edu",
  cqstools   => "/home/shengq1/cqstools/cqstools.exe",

  #mapping
  aligner   => "bwa",
  bwa_fasta => "/scratch/cqs/shengq1/references/gatk/b37/bwa_index_0.7.12/human_g1k_v37.fasta",

  #call variants
  perform_gatk_callvariants   => 1,
  gatk_callvariants_vqsr_mode => 0,
  gatk_jar                    => "/scratch/cqs/shengq1/local/bin/gatk/GenomeAnalysisTK.jar",
  picard_jar                  => "/scratch/cqs/shengq1/local/bin/picard/picard.jar",
  dbsnp                       => "/scratch/cqs/shengq1/references/gatk/b37/dbsnp_150.b37.vcf",
  hapmap                      => "/scratch/cqs/shengq1/references/gatk/b37/hapmap_3.3.b37.vcf",
  omni                        => "/scratch/cqs/shengq1/references/gatk/b37/1000G_omni2.5.b37.vcf",
  g1000                       => "/scratch/cqs/shengq1/references/gatk/b37/1000G_phase1.snps.high_confidence.b37.vcf",
  mills                       => "/scratch/cqs/shengq1/references/gatk/b37/Mills_and_1000G_gold_standard.indels.b37.vcf",

  #somatic mutation
  perform_muTect      => 1,
  muTect_init_command => "setpkgs -a java",
  muTect_option       => "--min_qscore 20 --filter_reads_with_N_cigar",
  muTect_jar          => "/home/shengq1/local/bin/mutect-1.1.7.jar",

  #annotation
  annovar_param    => "-protocol refGene,avsnp147,cosmic70,exac03 -operation g,f,f,f",
  annovar_db       => "/scratch/cqs/shengq1/references/annovar/humandb/",
  annovar_buildver => "hg19",

  files => {
    "TP0097_Norecur_Base1" => [
      "/scratch/h_vangard_1/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233961/BHCHJNBBXX_DS-233961_AGTTCC_L005_R1_001.fastq.gz",
      "/scratch/h_vangard_1/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233961/BHCHJNBBXX_DS-233961_AGTTCC_L005_R2_001.fastq.gz"
    ],
    "TP0109_Recur_Base1" => [
      "/scratch/h_vangard_1/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233949/AH7YCFBBXX_DS-233949_CGATGT_L005_R1_001.fastq.gz",
      "/scratch/h_vangard_1/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233949/AH7YCFBBXX_DS-233949_CGATGT_L005_R2_001.fastq.gz",
    ],
  },

  groups => {
    "TP0097_TP0109" => [ "TP0097_Norecur_Base1", "TP0109_Recur_Base1" ]
  }
};

performExomeSeq($def);

1;

