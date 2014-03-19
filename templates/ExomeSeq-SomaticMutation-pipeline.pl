#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ClassFactory;

my $vangard = "template";

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/temp/exome_somatic_mutation");

my $bwa_dir = "${target_dir}/bwa_refine";

my $fasta_file           = "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa";
my $transcript_gtf       = "/scratch/cqs/shengq1/references/hg19/Homo_sapiens.GRCh37.73.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_73";

#don't export as csv format since cqstools will only compatible with tab delimited format
my $annovar_param = "-protocol refGene,snp137,cosmic64,esp6500si_all,1000g2012apr_all -operation g,f,f,f,f --remove";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";

my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => "${vangard}" },
  fastqfiles => {
    "2055-PM-00" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-0_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-0_2.fastq.gz" ],
    "2055-PM-01" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-1_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-1_2.fastq.gz" ],
    "2055-PM-02" => [ "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-2_1.fastq.gz",  "/autofs/blue_sequencer/Runs/projects/2055-PM/2013-09-24/2055-PM-2_2.fastq.gz" ],
  },
  groups => { #should be normal, tumor pair
    "P4413_1_NL"       => [ "2055-PM-00", "2055-PM-01" ],
    "P4413_2_SEV_D"    => [ "2055-PM-00", "2055-PM-02" ],
  },
  fastqc => {
    class      => "FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqfiles",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa => {
    class      => "BWA",
    perform    => 1,
    target_dir => "${target_dir}/bwa",
    option     => "-q 15 -t 8",
    fasta_file => $fasta_file,
    source_ref => "fastqfiles",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  refine => {
    class              => "GATKRefine",
    perform            => 1,
    target_dir         => "${target_dir}/refine",
    option             => "-Xmx40g",
    fasta_file         => $fasta_file,
    source_ref         => "bwa",
    thread_count       => 8,
    vcf_files          => ["/data/cqs/shengq1/reference/snp137/human/00-All.vcf"],
    gatk_jar           => "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar",
    markDuplicates_jar => "/home/shengq1/local/bin/picard/MarkDuplicates.jar",
    sh_direct          => 0,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  muTect => {
    class       => "MuTect",
    perform     => 1,
    target_dir  => "${target_dir}/muTect",
    option      => "", #don't use thread mode, it may cause dead-lock
    source_ref  => "refine",
    groups_ref  => "groups",
    java_option => "-Xmx40g",
    fasta_file  => $fasta_file,
    cosmic_file => "/data/cqs/shengq1/reference/cosmic/cosmic_v65_28052013.hg19.16571.vcf",
    dbsnp_file  => "/data/cqs/shengq1/reference/snp137/human/00-All.vcf",
    sh_direct   => 1,
    muTect_jar  => "/home/shengq1/local/bin/muTect-1.1.4.jar",
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  annovar_mutect => {
    class      => "Annovar",
    perform    => 1,
    target_dir => "${target_dir}/muTect",
    option     => $annovar_param,
    source_ref => [ "muTect", ".pass.vcf\$" ],
    annovar_db => $annovar_db,
    buildver   => "hg19",
    cqstools   => $cqstools, #use cqstools to generate final excel report file
    affy_file  => "/data/cqs/shengq1/reference/affy/HG-U133_Plus_2.na33.annot.csv", #use affy file to get gene description information
    sh_direct  => 0,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);

1;
