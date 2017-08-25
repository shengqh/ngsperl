#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use Pipeline::RNASeq;

my $def = {
  task_name      => "B3436",
  cqstools       => "/home/shengq1/cqstools/cqstools.exe",
  email          => "quanhu.sheng\@vanderbilt.edu",
  target_dir     => create_directory_or_die("/scratch/cqs/shengq1/test/pipeline_rnaseq_simple"),
  transcript_gtf => "/scratch/cqs/shengq1/references/ensembl/v78/Mus_musculus.GRCm38.78.MT.gtf",
  name_map_file  => "/scratch/cqs/shengq1/references/ensembl/v78/Mus_musculus.GRCm38.78.MT.map",
  star_index     => "/scratch/cqs/shengq1/references/mm10_sorted_M/STAR_index_v38.78_2.5.0b_sjdb49",
  fasta_file     => "/scratch/cqs/shengq1/references/mm10_sorted_M/mm10.fa",
  files          => {
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

  perform_rnaseqc => 1,
  rnaseqc_jar     => "/home/shengq1/local/bin/RNA-SeQC_v1.1.8.jar",

};

performRNASeq($def);
1;

