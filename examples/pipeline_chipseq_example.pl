#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;
use Pipeline::SmallRNAUtils;
use Pipeline::ChIPSeq;

my $def = {
  task_name  => "Tessa-Ramos-ChIPseq",
  cqstools   => "/home/shengq2/cqstools/cqstools.exe",
  email      => "jing.wang.1\@vumc.org",
  target_dir => create_directory_or_die("/scratch/h_vangard_1/wangj52/Bill/Tessa-Ramos-ChIPseq"),
  picard_jar => "/scratch/cqs/shengq2/local/bin/picard/picard.jar",
  add_folder_index => 1,
  constraint => "haswell",
  pairend => 1,

  perform_cutadapt => 0,
  adapter          => "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",    #trueseq adapter
  min_read_length  => 30,
  
#  perform_cutadapt => 1,
#  cutadapt_option  => "-q 20 -a AGATCGGAAGAGCACACGTC -A AGATCGGAAGAGCGTCGTGT",
#  min_read_length  => 30,


  aligner       => "bowtie2",
#  bowtie1_fasta => "/scratch/cqs/shengq2/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome.fa",
  bowtie2_index => "/scratch/cqs/references/human/hg19/bowtie2_index_2.3.4.1/GRCh37.p13.genome",

  #peak
  #peak_caller => "macs1",

  peak_caller     => "macs2",
  macs2_genome    => "hs",
  macs2_peak_type => "narrow",
  macs2_option    => "-B -q 0.05 -g hs -f BAMPE",

  #macs2_peak_type => "broad",
  #macs2_option => "--broad --broad-cutoff 0.1 -B -q 0.01 -g hs",
  
  perform_homer_motifs => 1,
  homer_option => "",
  homer_genome => "hg19",
  
  
  perform_rose => 0,
  rose_folder  => "/scratch/cqs/shengq2/local/bin/bradnerlab",
  rose_genome  => "hg19",

  perform_bamplot => 0,
  bamplot_gff     => "/scratch/cqs/shengq2/brown/20170116_chipseq_3593_11to15_human/document/HAEC.gff",
  dataset_name    => "HG19",
  bamplot_option  => "-g HG19 -y uniform -r",

  files => {
    "MYC_HA_N1"       => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-1-ATCACG_S1_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-1-ATCACG_S1_R2_001.fastq.gz"],
    "MYC_HA_N2"       => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-2-ACAGTG_S2_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-2-ACAGTG_S2_R2_001.fastq.gz"],
    "MYC_HA_N3"       => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-3-GATCAG_S3_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-3-GATCAG_S3_R2_001.fastq.gz"],
    "MYC_WT_N1"       => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-4-CGATGT_S4_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-4-CGATGT_S4_R2_001.fastq.gz"],
    "MYC_WT_N2"       => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-5-GCCAAT_S5_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-5-GCCAAT_S5_R2_001.fastq.gz"],
    "MYC_WT_N3"       => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-6-TAGCTT_S6_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-6-TAGCTT_S6_R2_001.fastq.gz"],
    "MYC_4A_N1"       => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-7-TTAGGC_S1_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-7-TTAGGC_S1_R2_001.fastq.gz"],
    "MYC_4A_N2"       => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-8-CAGATC_S2_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-8-CAGATC_S2_R2_001.fastq.gz"],
    "MYC_4A_N3"       => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-9-GGCTAC_S3_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-9-GGCTAC_S3_R2_001.fastq.gz"],
    "MYC_EHAY_N1"       => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-10-TGACCA_S4_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-10-TGACCA_S4_R2_001.fastq.gz"],
    "MYC_EHAY_N2"       => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-11-ACTTGA_S5_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-11-ACTTGA_S5_R2_001.fastq.gz"],
    "MYC_EHAY_N3"       => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-12-CTTGTA_S6_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/2889/2889-TP-12-CTTGTA_S6_R2_001.fastq.gz"],
    "HCF1_IgG_N1"        => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-1-ATCACG_S1_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-1-ATCACG_S1_R2_001.fastq.gz"],
    "HCF1_IgG_N2"        => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-2-ACAGTG_S2_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-2-ACAGTG_S2_R2_001.fastq.gz"],
    "HCF1_IgG_N3"        => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-3-GATCAG_S3_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-3-GATCAG_S3_R2_001.fastq.gz"],
    "HCF1_WT_N1"        => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-4-CGATGT_S4_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-4-CGATGT_S4_R2_001.fastq.gz"],
    "HCF1_WT_N2"        => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-5-GCCAAT_S5_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-5-GCCAAT_S5_R2_001.fastq.gz"],
    "HCF1_WT_N3"        => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-6-ACTGAT_S6_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-6-ACTGAT_S6_R2_001.fastq.gz"],
    "HCF1_4A_N1"        => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-7-TTAGGC_S1_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-7-TTAGGC_S1_R2_001.fastq.gz"],
    "HCF1_4A_N2"        => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-8-CAGATC_S2_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-8-CAGATC_S2_R2_001.fastq.gz"],
    "HCF1_4A_N3"        => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-9-GGCTAC_S3_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-9-GGCTAC_S3_R2_001.fastq.gz"],
    "HCF1_EHAY_N1"        => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-10-TGACCA_S4_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-10-TGACCA_S4_R2_001.fastq.gz"],
    "HCF1_EHAY_N2"        => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-11-ACTTGA_S5_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-11-ACTTGA_S5_R2_001.fastq.gz"],
    "HCF1_EHAY_N3"        => ["/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-12-CTTGTA_S6_R1_001.fastq.gz", "/data/cqs/wangj52/Bill/Tessa-Ramos-ChIPseq/3107/3107-TP-12-CTTGTA_S6_R2_001.fastq.gz"],
  },
  treatments => {
    "MYC_WT_N1"         => ["MYC_WT_N1"],
    "MYC_WT_N2"         => ["MYC_WT_N2"],
    "MYC_WT_N3"         => ["MYC_WT_N3"],
    "MYC_4A_N1"         => ["MYC_4A_N1"],
    "MYC_4A_N2"         => ["MYC_4A_N2"],
    "MYC_4A_N3"         => ["MYC_4A_N3"],
    "MYC_EHAY_N1"         => ["MYC_EHAY_N1"],
    "MYC_EHAY_N2"         => ["MYC_EHAY_N2"],
    "MYC_EHAY_N3"         => ["MYC_EHAY_N3"],
    "HCF1_WT_N1"         => ["HCF1_WT_N1"],
    "HCF1_WT_N2"         => ["HCF1_WT_N2"],
    "HCF1_WT_N3"         => ["HCF1_WT_N3"],
    "HCF1_4A_N1"         => ["HCF1_4A_N1"],
    "HCF1_4A_N2"         => ["HCF1_4A_N2"],
    "HCF1_4A_N3"         => ["HCF1_4A_N3"],
    "HCF1_EHAY_N1"         => ["HCF1_EHAY_N1"],
    "HCF1_EHAY_N2"         => ["HCF1_EHAY_N2"],
    "HCF1_EHAY_N3"         => ["HCF1_EHAY_N3"],
  },
  controls => {
    "MYC_WT_N1"         => ["MYC_HA_N1"],
    "MYC_WT_N2"         => ["MYC_HA_N2"],
    "MYC_WT_N3"         => ["MYC_HA_N3"],
    "MYC_4A_N1"         => ["MYC_HA_N1"],
    "MYC_4A_N2"         => ["MYC_HA_N2"],
    "MYC_4A_N3"         => ["MYC_HA_N3"],
    "MYC_EHAY_N1"         => ["MYC_HA_N1"],
    "MYC_EHAY_N2"         => ["MYC_HA_N2"],
    "MYC_EHAY_N3"         => ["MYC_HA_N3"],
    "HCF1_WT_N1"         => ["HCF1_IgG_N1"],
    "HCF1_WT_N2"         => ["HCF1_IgG_N2"],
    "HCF1_WT_N3"         => ["HCF1_IgG_N3"],
    "HCF1_4A_N1"         => ["HCF1_IgG_N1"],
    "HCF1_4A_N2"         => ["HCF1_IgG_N2"],
    "HCF1_4A_N3"         => ["HCF1_IgG_N3"],
    "HCF1_EHAY_N1"         => ["HCF1_IgG_N1"],
    "HCF1_EHAY_N2"         => ["HCF1_IgG_N2"],
    "HCF1_EHAY_N3"         => ["HCF1_IgG_N3"],
  },

  perform_chipqc          => 1,
  chipqc_genome           => "hg19",
  perform_diffbind        => 1,
  homer_annotation_genome => "hg19",
  chipqc_chromosomes      => 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22',
  design_table            => {
    Tissue      => "Ramos",
    "MYC" => {
    	Factor    => "MYC", 
      "MYC_WT_N1" => {
        Condition => "MYC_WT",
        Replicate => "1"
      },
      "MYC_WT_N2" => {
        Condition => "MYC_WT",
        Replicate => "2"
      },
      "MYC_WT_N3" => {
        Condition => "MYC_WT",
        Replicate => "3"
      },
      "MYC_4A_N1" => {
        Condition => "MYC_4A",
        Replicate => "1"
      },
      "MYC_4A_N2" => {
        Condition => "MYC_4A",
        Replicate => "2"
      },
      "MYC_4A_N3" => {
        Condition => "MYC_4A",
        Replicate => "3"
      },
      "MYC_EHAY_N1" => {
        Condition => "MYC_EHAY",
        Replicate => "1"
      },
      "MYC_EHAY_N2" => {
        Condition => "MYC_EHAY",
        Replicate => "2"
      },
      "MYC_EHAY_N3" => {
        Condition => "MYC_EHAY",
        Replicate => "3"
      },
      "Comparison" => [
         [ "MYC-4A_vs_WT", "MYC_WT", "MYC_4A" ],
         [ "MYC-EHAY_vs_WT", "MYC_WT", "MYC_EHAY" ],
       ],
       "MinOverlap" => {
         "MYC_WT" => 2,
         "MYC_4A" => 2,
         "MYC_EHAY" => 2,
      },
    },
    "HCF_1" => {
    	Factor    => "HCF_1", 
            "HCF1_WT_N1" => {
        Condition => "HCF1_WT",
        Replicate => "1"
      },
      "HCF1_WT_N2" => {
        Condition => "HCF1_WT",
        Replicate => "2"
      },
      "HCF1_WT_N3" => {
        Condition => "HCF1_WT",
        Replicate => "3"
      },
      "HCF1_4A_N1" => {
        Condition => "HCF1_4A",
        Replicate => "1"
      },
      "HCF1_4A_N2" => {
        Condition => "HCF1_4A",
        Replicate => "2"
      },
      "HCF1_4A_N3" => {
        Condition => "HCF1_4A",
        Replicate => "3"
      },
      "HCF1_EHAY_N1" => {
        Condition => "HCF1_EHAY",
        Replicate => "1"
      },
      "HCF1_EHAY_N2" => {
        Condition => "HCF1_EHAY",
        Replicate => "2"
      },
      "HCF1_EHAY_N3" => {
        Condition => "HCF1_EHAY",
        Replicate => "3"
      },
      "Comparison" => [
         [ "HCF1-4A_vs_WT", "HCF1_WT", "HCF1_4A" ],
         [ "HCF1-EHAY_vs_WT", "HCF1_WT", "HCF1_EHAY" ],
       ],
       "MinOverlap" => {
         "HCF1_WT" => 2,
         "HCF1_4A" => 2,
         "HCF1_EHAY" => 2,
      },
    },
  },
};

performChIPSeq($def);
1;
