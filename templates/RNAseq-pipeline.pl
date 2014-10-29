#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/pipelines/RNAseq-pipeline");

my $fasta_file           = "/data/cqs/guoy1/reference/mm10/bowtie2_index/mm10.fa";
my $bowtie2_index        = "/data/cqs/guoy1/reference/mm10/bowtie2_index/mm10";
my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/mm10/Mus_musculus.GRCm38.68_chr1-22-X-Y-M.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/references/mm10/gtfindex/Mus_musculus.GRCm38.68";
my $dexseq_gff           = "/data/cqs/guoy1/reference/annotation2/mm10/Mus_musculus.GRCm38.74_chr1-19-X-Y-M.dexseq.gff";
my $dexseqpy             = "/home/shengq1/pylibs/bin/dexseq_count.py";
my $name_map_file        = "/data/cqs/guoy1/reference/annotation2/mm10/Mus_musculus.GRCm38.74_chr1-19-X-Y-M.map";
my $files                = {
  "S1" => ["/gpfs21/scratch/cqs/shengq1/report/rawdata/s1_sequence.txt"],
  "S2" => ["/gpfs21/scratch/cqs/shengq1/report/rawdata/s2_sequence.txt"],
  "S3" => ["/gpfs21/scratch/cqs/shengq1/report/rawdata/s3_sequence.txt"],
  "S4" => ["/gpfs21/scratch/cqs/shengq1/report/rawdata/s4_sequence.txt"],
  "S5" => ["/gpfs21/scratch/cqs/shengq1/report/rawdata/s5_sequence.txt"],
  "S6" => ["/gpfs21/scratch/cqs/shengq1/report/rawdata/s6_sequence.txt"],
};
my $groups = {
  "G1" => [ "S1", "S2", "S3" ],
  "G2" => [ "S4", "S5", "S6" ],
};
my $pairs    = { "G2_vs_G1" => [ "G1", "G2" ] };
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";
my $email    = "quanhu.sheng\@vanderbilt.edu";
my $task     = "pipeline";

my $def = {
  tophat_option => "--segment-length 25 -r 0 -p 6",
};

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
    option               => $def->{tophat_option},
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
  rnaseqc => {
    class          => "RNASeQC",
    perform        => 1,
    target_dir     => "${target_dir}/RNASeQC",
    option         => "",
    fasta_file     => $fasta_file,
    transcript_gtf => $transcript_gtf,
    jar            => "/home/shengq1/local/bin/RNA-SeQC_v1.1.7.jar",
    source_ref     => "tophat2",
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  cufflinks => {
    class          => "Cufflinks",
    perform        => 1,
    target_dir     => "${target_dir}/cufflinks",
    option         => "-p 8 -u -N",
    transcript_gtf => $transcript_gtf,
    source_ref     => "tophat2",
    sh_direct      => 1,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  cuffmerge => {
    class      => "Cuffmerge",
    perform    => 1,
    target_dir => "${target_dir}/cuffmerge",
    option     => "-p 8",
    source_ref => "cufflinks",
    fasta_file => $fasta_file,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  cufflinks_cuffdiff => {
    class              => "Cuffdiff",
    perform            => 1,
    target_dir         => "${target_dir}/cufflinks_cuffdiff",
    option             => "-p 8 -u -N",
    transcript_gtf_ref => "cuffmerge",
    source_ref         => "tophat2",
    groups             => $groups,
    pairs              => $pairs,
    fasta_file         => $fasta_file,
    sh_direct          => 1,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "720",
      "mem"      => "40gb"
    },
  },
  cuffdiff => {
    class          => "Cuffdiff",
    perform        => 1,
    target_dir     => "${target_dir}/cuffdiff",
    option         => "-p 8 -u -N",
    transcript_gtf => $transcript_gtf,
    source_ref     => "tophat2",
    groups         => $groups,
    pairs          => $pairs,
    fasta_file     => $fasta_file,
    sh_direct      => 1,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "720",
      "mem"      => "40gb"
    },
  },
  sortbam => {
    class         => "Sortbam",
    perform       => 1,
    target_dir    => "${target_dir}/sortname",
    option        => "",
    source_ref    => "tophat2",
    sort_by_query => 1,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  htseqcount => {
    class      => "HTSeqCount",
    perform    => 1,
    target_dir => "${target_dir}/htseqcount",
    option     => "",
    source_ref => "sortbam",
    gff_file   => $transcript_gtf,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  genetable => {
    class         => "CQSDatatable",
    perform       => 1,
    target_dir    => "${target_dir}/genetable",
    option        => "-p ENS --noheader -o ${task}_gene.count",
    source_ref    => "htseqcount",
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
  deseq2 => {
    class         => "DESeq2",
    perform       => 1,
    target_dir    => "${target_dir}/deseq2",
    option        => "",
    source        => $pairs,
    groups        => $groups,
    countfile_ref => "genetable",
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  dexseqcount => {
    class        => "DexseqCount",
    perform      => 1,
    target_dir   => "${target_dir}/dexseqcount",
    option       => "",
    source_ref   => "tophat2",
    gff_file     => $transcript_gtf,
    dexseq_count => $dexseqpy,
    sh_direct    => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  exontable => {
    class         => "CQSDatatable",
    perform       => 1,
    target_dir    => "${target_dir}/exontable",
    option        => "-p ENS --noheader -o ${task}_exon.count",
    name_map_file => $name_map_file,
    source_ref    => "dexseqcount",
    cqs_tools     => $cqstools,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  sequencetask => {
    class      => "SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    source     => {
      "sample" => [ "fastqc",  "tophat2",   "cufflinks",          "sortbam",  "htseqcount", "dexseqcount" ],
      "task"   => [ "rnaseqc", "cuffmerge", "cufflinks_cuffdiff", "cuffdiff", "genetable",  "deseq2", "exontable" ],
    },
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "30gb"
    }
  }
};

performConfig($config);

1;
