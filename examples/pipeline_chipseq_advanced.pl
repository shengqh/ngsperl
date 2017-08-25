#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/chipseq/pipeline");

my $fasta_file     = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome.fa";
my $bowtie_index   = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome";
my $cqstools       = "/home/shengq1/cqstools/cqstools.exe";
my $qc3_perl       = "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl";
my $transcript_gtf = "/scratch/cqs/shengq1/references/gencode/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "195R3";

my $config = {
  general => { task_name => $task },
  files   => {
    "d17_static_ESctrl2_Input"   => ["/gpfs21/scratch/cqs/shengq1/chipseq/20160302_janathan_chipseq_195R3/195R/day_17_static/CT663p1-4_bc01_ATCACG_L006_R1_001.fastq.gz"],
    "d17_static_ESctrl2_H3K27ac" => ["/gpfs21/scratch/cqs/shengq1/chipseq/20160302_janathan_chipseq_195R3/195R/day_17_static/CT663p55-62_index11_GGCTAC_L005_R1_001.fastq.gz"],
    "d17_shear_ESctrl2_Input"    => ["/gpfs21/scratch/cqs/shengq1/chipseq/20160302_janathan_chipseq_195R3/195R/day_17_shear/CT663p1-4_bc03_TTAGGC_L006_R1_001.fastq.gz"],
    "d17_shear_ESctrl2_H3K27ac"  => ["/gpfs21/scratch/cqs/shengq1/chipseq/20160302_janathan_chipseq_195R3/195R/day_17_shear/CT663p55-62_index7_CAGATC_L005_R1_001.fastq.gz"],
  },
  treatments => {
    "static" => ["d17_static_ESctrl2_H3K27ac"],
    "shear"  => ["d17_shear_ESctrl2_H3K27ac"],
  },
  controls => {
    "static" => ["d17_static_ESctrl2_Input"],
    "shear"  => ["d17_shear_ESctrl2_Input"],
  },
  diffpairs => {
    "shear_vs_static" => [ "static", "shear" ],
  },
  fastqc_pre_trim => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_pre_trim",
    option     => "",
    source_ref => ["files"],
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc_pre_trim_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    sh_direct  => 1,
    target_dir => "${target_dir}/fastqc_pre_trim",
    cqstools   => $cqstools,
    option     => "",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  cutadapt => {
    class      => "Trimmer::Cutadapt",
    perform    => 1,
    target_dir => "${target_dir}/cutadapt",
    option     => "-m 30 --trim-n",
    source_ref => ["files"],
    adapter    => "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",    #trueseq adapter
    extension  => "_clipped.fastq",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  fastqc_post_trim => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc_post_trim",
    option     => "",
    sh_direct  => 1,
    source_ref => [ "cutadapt", ".fastq.gz" ],
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc_post_trim_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    sh_direct  => 1,
    target_dir => "${target_dir}/fastqc_post_trim",
    cqstools   => $cqstools,
    option     => "",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastq_len => {
    class      => "CQS::FastqLen",
    perform    => 1,
    target_dir => "$target_dir/fastq_len",
    option     => "",
    source_ref => "cutadapt",
    cqstools   => $cqstools,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bowtie1 => {
    class                   => "Alignment::Bowtie1",
    perform                 => 1,
    target_dir              => "${target_dir}/bowtie1",
    option                  => "-v 1 -m 1 --best --strata",
    fasta_file              => $fasta_file,
    source_ref              => [ "cutadapt", ".fastq.gz\$" ],
    bowtie1_index           => $bowtie_index,
    chromosome_grep_pattern => "\"^chr\"",
    sh_direct               => 0,
    pbs                     => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs1callpeak => {
    class        => "Chipseq::MACS",
    perform      => 1,
    target_dir   => "${target_dir}/macs1callpeak",
    option       => "-p 1e-9 -w -S --space=50",
    source_ref   => "bowtie1",
    groups_ref   => "treatments",
    controls_ref => "controls",
    sh_direct    => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs1callpeak_bradner_rose2 => {
    class                => "Chipseq::BradnerRose2",
    perform              => 1,
    target_dir           => "${target_dir}/macs1callpeak_bradner_rose2",
    option               => "",
    source_ref           => "bowtie1",
    groups_ref           => "treatments",
    controls_ref         => "controls",
    pipeline_dir         => "/scratch/cqs/shengq1/local/bin/bradnerlab",
    binding_site_bed_ref => [ "macs1callpeak", ".bed\$" ],
    genome               => "hg19",
    sh_direct            => 1,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs2callpeak => {
    class        => "Chipseq::MACS2Callpeak",
    perform      => 1,
    target_dir   => "${target_dir}/macs2callpeak",
    option       => "-g hs --broad -B -p 1e-9",
    source_ref   => "bowtie1",
    groups_ref   => "treatments",
    controls_ref => "controls",
    sh_direct    => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs2bdgdiff => {
    class      => "Chipseq::MACS2Bdgdiff",
    perform    => 1,
    target_dir => "${target_dir}/macs2bdgdiff",
    option     => "",
    source_ref => "macs2callpeak",
    groups_ref => "diffpairs",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_1 => [ "fastqc_pre_trim",         "cutadapt",                 "fastqc_post_trim", "bowtie1",                     "fastq_len" ],
      step_2 => [ "fastqc_pre_trim_summary", "fastqc_post_trim_summary", "macs1callpeak",    "macs1callpeak_bradner_rose2", "macs2callpeak" ],
      step_3 => ["macs2bdgdiff"],
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
