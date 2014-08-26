#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/pipelines/RNAseq-pipeline-pairend");

my $fasta_file           = "/data/cqs/guoy1/reference/mm10/bowtie2_index/mm10.fa";
my $bowtie2_index        = "/data/cqs/guoy1/reference/mm10/bowtie2_index/mm10";
my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/mm10/Mus_musculus.GRCm38.68_chr1-22-X-Y-M.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/references/mm10/gtfindex/Mus_musculus.GRCm38.68";
my $dexseq_gff           = "/data/cqs/guoy1/reference/annotation2/mm10/Mus_musculus.GRCm38.74_chr1-19-X-Y-M.dexseq.gff";
my $dexseqpy             = "/home/shengq1/pylibs/bin/dexseq_count.py";
my $name_map_file        = "/data/cqs/guoy1/reference/annotation2/mm10/Mus_musculus.GRCm38.74_chr1-19-X-Y-M.map";
my $files                = {
  "S12" => ["/gpfs21/scratch/cqs/shengq1/report/rawdata/s1_sequence.txt,/gpfs21/scratch/cqs/shengq1/report/rawdata/s2_sequence.txt"],
  "S34" => ["/gpfs21/scratch/cqs/shengq1/report/rawdata/s3_sequence.txt,/gpfs21/scratch/cqs/shengq1/report/rawdata/s4_sequence.txt"],
  "S56" => ["/gpfs21/scratch/cqs/shengq1/report/rawdata/s5_sequence.txt,/gpfs21/scratch/cqs/shengq1/report/rawdata/s6_sequence.txt"],
};
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";
my $email    = "quanhu.sheng\@vanderbilt.edu";
my $task     = "pipeline";

my $config = {
  general => { task_name => $task },
  fastqc  => {
    class      => "FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source     => $files,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  tophat2 => {
    class                => "Tophat2",
    perform              => 1,
    target_dir           => "${target_dir}/tophat2",
    option               => "--segment-length 25 -r 0 -p 6",
    source               => $files,
    bowtie2_index        => $bowtie2_index,
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    rename_bam           => 1,
    sh_direct            => 1,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=6",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
};

performConfig($config);

1;
