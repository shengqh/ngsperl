#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;
use Pipeline::SmallRNAUtils;
use Pipeline::ChIPSeq;

my $def = {
  task_name  => "zf_264_1to2_mouse",
  cqstools   => "/home/shengq1/cqstools/cqstools.exe",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => create_directory_or_die("/scratch/cqs/shengq1/brown/20170530_chipseq_zf_264_1to2_mouse"),

  files => {
    "H3K27ac" => ["/scratch/cqs/shengq1/brown/data/264/264-ZF-1_S1_R1_001.fastq.gz"],
    "Input"   => ["/scratch/cqs/shengq1/brown/data/264/264-ZF-2_S2_R1_001.fastq.gz"],
  },
  treatments => {
    "H3K27ac" => ["H3K27ac"],
  },
  controls => {
    "H3K27ac" => ["Input"],
  },
  
  add_folder_index => 1,

  #preprocessing
  perform_cutadapt => 1,
  adapter          => "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",    #trueseq adapter
  min_read_length  => 30,

  #mapping
  aligner       => "bowtie2",
  bowtie2_index => "/scratch/cqs/shengq1/references/ucsc/illumina/mm10/Sequence/Bowtie2Index/genome",
  picard_jar    => "/scratch/cqs/shengq1/local/bin/picard/picard.jar",

  #peak calling
  peak_caller => "macs",

  #homer
  perform_homer_motifs => 1,
  homer_option => "",
  homer_genome => "mm10",

  #annotation
  perform_rose => 0,
  rose_folder  => "/scratch/cqs/shengq1/local/bin/bradnerlab",
  rose_genome  => "mm10",

  perform_bamplot => 0,
  bamplot_gff     => "/scratch/cqs/shengq1/brown/20170116_chipseq_3593_11to15_human/document/HAEC.gff",
  dataset_name    => "MM10",
  bamplot_option  => "-g MM10 -y uniform -r",

  #qc
  perform_chipqc     => 0,
  chipqc_genome      => "mm10",
  chipqc_chromosomes => 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19',
  design_table       => {
    "H3K27ac" => {
      "H3K27ac" => {
        Condition => "H3K27ac",
        Replicate => "1"
      },
    },
  },

};

performChIPSeq($def);
1;
