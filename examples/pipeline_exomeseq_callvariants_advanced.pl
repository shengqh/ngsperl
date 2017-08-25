#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/test/20161013_liuqi_gene_panel");
my $cqstools   = "/home/shengq1/cqstools/cqstools.exe";
my $email      = "quanhu.sheng.1\@vanderbilt.edu";

my $bwa_fasta = "/scratch/cqs/shengq1/references/gatk/b37/bwa_index_0.7.12/human_g1k_v37.fasta";

#my $dbsnp  = "/scratch/cqs/shengq1/references/dbsnp/human_9606_b150_GRCh37p13.vcf";
my $dbsnp  = "/scratch/cqs/shengq1/references/gatk/b37/dbsnp_150.b37.vcf";
my $hapmap = "/scratch/cqs/shengq1/references/gatk/b37/hapmap_3.3.b37.vcf";
my $omni   = "/scratch/cqs/shengq1/references/gatk/b37/1000G_omni2.5.b37.vcf";
my $g1000  = "/scratch/cqs/shengq1/references/gatk/b37/1000G_phase1.snps.high_confidence.b37.vcf";
my $mills  = "/scratch/cqs/shengq1/references/gatk/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";

my $cosmic = "/scratch/cqs/shengq1/references/cosmic/cosmic_v71_hg19_16569_MT.vcf";

my $annovar_param = "-protocol refGene,avsnp147,cosmic70,exac03 -operation g,f,f,f";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";
my $gatk_jar      = "/scratch/cqs/shengq1/local/bin/gatk/GenomeAnalysisTK.jar";
my $picard_jar    = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";

my $task_name = "Adenoma";
my $config    = {
  general => { task_name => $task_name },
  files   => {
    "TP0097_Norecur_Base1" => [
      "/scratch/h_vangard_1/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233961/BHCHJNBBXX_DS-233961_AGTTCC_L005_R1_001.fastq.gz",
      "/scratch/h_vangard_1/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233961/BHCHJNBBXX_DS-233961_AGTTCC_L005_R2_001.fastq.gz"
    ],
    "TP0109_Recur_Base1" => [
      "/scratch/h_vangard_1/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233949/AH7YCFBBXX_DS-233949_CGATGT_L005_R1_001.fastq.gz",
      "/scratch/h_vangard_1/liuq6/Bob/Adenoma/RAW_DATA/Sample_DS-233949/AH7YCFBBXX_DS-233949_CGATGT_L005_R2_001.fastq.gz",
    ],
  },
  fastqc => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "files",
    sh_direct  => 1,
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
    cqstools   => $cqstools,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa => {
    class      => "Alignment::BWA",
    perform    => 1,
    target_dir => "${target_dir}/bwa",
    option     => "",
    bwa_index  => $bwa_fasta,
    source_ref => "files",
    picard_jar => $picard_jar,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  },
  bwa_refine => {
    class      => "GATK::Refine",
    perform    => 1,
    target_dir => "${target_dir}/bwa_refine",
    option     => "-Xmx40g",

    #gatk_option => "--fix_misencoded_quality_scores",
    gatk_option              => "",
    fasta_file               => $bwa_fasta,
    source_ref               => "bwa",
    vcf_files                => [ $dbsnp, $mills ],
    gatk_jar                 => $gatk_jar,
    picard_jar               => $picard_jar,
    sh_direct                => 0,
    slim_print_reads         => 1,
    samtools_baq_calibration => 0,
    sorted                   => 1,
    pbs                      => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf => {
    class         => "GATK::HaplotypeCaller",
    perform       => 1,
    target_dir    => "${target_dir}/bwa_refine_hc_gvcf",
    option        => "",
    source_ref    => "bwa_refine",
    java_option   => "",
    fasta_file    => $bwa_fasta,
    gatk_jar      => $gatk_jar,
    extension     => ".g.vcf",
    by_chromosome => 0,
    gvcf          => 1,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf_vqsr => {
    class       => "GATK::VariantFilter",
    perform     => 1,
    target_dir  => "${target_dir}/bwa_refine_hc_gvcf_vqsr",
    option      => "",
    vqsr_mode   => 1,
    source_ref  => "bwa_refine_hc_gvcf",
    java_option => "",
    fasta_file  => $bwa_fasta,
    dbsnp_vcf   => $dbsnp,
    hapmap_vcf  => $hapmap,
    omni_vcf    => $omni,
    g1000_vcf   => $g1000,
    mills_vcf   => $mills,
    gatk_jar    => $gatk_jar,
    cqstools    => $cqstools,
    sh_direct   => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf_vqsr_annovar => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/bwa_refine_hc_gvcf_vqsr_annovar",
    source_ref => "bwa_refine_hc_gvcf_vqsr",
    option     => $annovar_param,
    annovar_db => $annovar_db,
    buildver   => "hg19",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa_refine_hc_gvcf_hardfilter => {
    class       => "GATK::VariantFilter",
    perform     => 1,
    target_dir  => "${target_dir}/bwa_refine_hc_gvcf_hardfilter",
    option      => "",
    source_ref  => "bwa_refine_hc_gvcf",
    java_option => "",
    gatk_jar    => $gatk_jar,
    fasta_file  => $bwa_fasta,
    sh_direct   => 1,
    vqsr_mode   => 0,
    is_rna      => 0,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf_hardfilter_annovar => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/bwa_refine_hc_gvcf_hardfilter_annovar",
    source_ref => "bwa_refine_hc_gvcf_hardfilter",
    option     => $annovar_param,
    annovar_db => $annovar_db,
    buildver   => "hg19",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => [ "fastqc",         "bwa",                     "bwa_refine",                      "bwa_refine_hc_gvcf" ],
      step2 => [ "fastqc_summary", "bwa_refine_hc_gvcf_vqsr", "bwa_refine_hc_gvcf_vqsr_annovar", "bwa_refine_hc_gvcf_hardfilter", "bwa_refine_hc_gvcf_hardfilter_annovar" ],
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

