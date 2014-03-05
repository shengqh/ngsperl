#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root     = "/scratch/cqs/shengq1/mirna/20140305_pipeline";
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $hg19_mrna_gff     = "/data/cqs/shengq1/reference/miRBase20/hsa.gff3";
my $hg19_trna_bed     = "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed";
my $hg19_trna_fasta   = "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed.fa";
my $hg19_smallrna_bed = "/data/cqs/guoy1/reference/smallrna/hg19_smallRNA_ucsc_ensembl.bed";

my $target_dir = create_directory_or_die($root);

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "mirna";

my $samtools = "/home/shengq1/local/bin/samtools/samtools";

my $bowtie1_option_1mm = "-a -m 100 --best --strata -v 1 -l 12 -p 8";

my $bowtie1_human_index = "/data/cqs/guoy1/reference/hg19/bowtie_index/hg19";

my $mirnacount_option          = "-s";                                                    #ignore score
my $trnacount_option           = "--length --sequence";
my $mirna_overlap_count_option = "-s --gtf_key miRNA";
my $fasta_file                 = "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa";

my $def = {
  source => {
    "2245-CRF-00" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-0_1.fastq.gz"],
    "2245-CRF-01" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-1_1.fastq.gz"],
    "2245-CRF-02" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-2_1.fastq.gz"],
    "2245-CRF-03" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-3_1.fastq.gz"],
    "2245-CRF-04" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-4_1.fastq.gz"],
    "2245-CRF-05" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-5_1.fastq.gz"],
    "2245-CRF-06" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-6_1.fastq.gz"],
    "2245-CRF-07" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-7_1.fastq.gz"],
    "2245-CRF-08" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-8_1.fastq.gz"],
    "2245-CRF-09" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-9_1.fastq.gz"],
    "2245-CRF-10" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-10_1.fastq.gz"],
    "2245-CRF-11" => ["/autofs/blue_sequencer/Runs/projects/2245-CRF/2012-08-03/2245-CRF-11_1.fastq.gz"],
  },
  coordinate          => $hg19_mrna_gff,
  trna_coordinate     => $hg19_trna_bed,
  trna_fasta          => $hg19_trna_fasta,
  smallrna_coordinate => $hg19_smallrna_bed,
  bowtie1_index       => $bowtie1_human_index,
  target_dir          => $target_dir . "/human2245",
  task_name           => $task_name . "_human2245",
  groups              => {
    "2245-CRF-NORMAL"          => [ "2245-CRF-00", "2245-CRF-01", "2245-CRF-02", "2245-CRF-03" ],
    "2245-CRF-STEATOSIS"       => [ "2245-CRF-04", "2245-CRF-05", "2245-CRF-06", "2245-CRF-07" ],
    "2245-CRF-STEATOPEPATITIS" => [ "2245-CRF-08", "2245-CRF-09", "2245-CRF-10", "2245-CRF-11" ],
  },
};

my $cur_target_dir = create_directory_or_die( $def->{target_dir} );
my $config         = {
  general  => { "task_name" => $def->{task_name}, },
  cutadapt => {
    class      => "Cutadapt",
    perform    => 1,
    target_dir => "${cur_target_dir}/cutadapt",
    option     => "-O 10 -m 12",
    source     => $def->{source},
    adaptor    => "TGGAATTCTCGGGTGCCAAGG",
    extension  => "_clipped.fastq",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  fastqlen => {
    class      => "FastqLen",
    perform    => 1,
    target_dir => "${cur_target_dir}/fastqlen",
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
  identical => {
    class      => "FastqIdentical",
    perform    => 1,
    target_dir => "${cur_target_dir}/identical",
    option     => "",
    source_ref => [ "cutadapt", ".fastq.gz" ],
    cqstools   => $cqstools,
    extension  => "_clipped_identical.fastq",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },

  #not identical, for IGV
  bowtie1_genome_cutadapt_topN_1mm_notidentical => {
    class         => "Bowtie1",
    perform       => 1,
    target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_notidentical",
    option        => $bowtie1_option_1mm,
    source_ref    => [ "cutadapt", ".fastq.gz" ],
    bowtie1_index => $def->{bowtie1_index},
    samonly       => 0,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },

  #1 mismatch search
  bowtie1_genome_cutadapt_topN_1mm => {
    class         => "Bowtie1",
    perform       => 1,
    target_dir    => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm",
    option        => $bowtie1_option_1mm,
    source_ref    => [ "identical", ".fastq\$" ],
    bowtie1_index => $def->{bowtie1_index},
    samonly       => 0,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  mirna_1mm_count => {
    class           => "MirnaCount",
    perform         => 1,
    target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA",
    option          => $mirnacount_option,
    source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
    fastq_files_ref => "identical",
    seqcount_ref    => [ "identical", ".dupcount\$" ],
    cqs_tools       => $cqstools,
    gff_file        => $def->{coordinate},
    fasta_file      => $fasta_file,
    samtools        => $samtools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  miRNA_1mm_table => {
    class      => "CQSMirnaTable",
    perform    => 1,
    target_dir => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_table",
    option     => "",
    source_ref => "mirna_1mm_count",
    cqs_tools  => $cqstools,
    groups     => $def->{groups},
    prefix     => "miRNA_1mm_",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  miRNA_1mm_count_overlap => {
    class           => "CQSMappedCount",
    perform         => 1,
    target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap",
    option          => $mirna_overlap_count_option,
    source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
    fastq_files_ref => "identical",
    seqcount_ref    => [ "identical", ".dupcount\$" ],
    cqs_tools       => $cqstools,
    gff_file        => $def->{coordinate},
    fasta_file      => $fasta_file,
    samtools        => $samtools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  miRNA_1mm_overlap_position => {
    class      => "CQSMappedPosition",
    perform    => 1,
    target_dir => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_overlap_position",
    option     => "-o " . $def->{task_name} . "_miRNA.position",
    source_ref => "miRNA_1mm_count_overlap",
    cqs_tools  => $cqstools,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  tRNA_1mm_count => {
    class           => "CQSMappedCount",
    perform         => 1,
    target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA",
    option          => $trnacount_option,
    source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
    fastq_files_ref => "identical",
    seqcount_ref    => [ "identical", ".dupcount\$" ],
    cqs_tools       => $cqstools,
    gff_file        => $def->{trna_coordinate},
    fasta_file      => $def->{trna_fasta},
    samtools        => $samtools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  tRNA_1mm_table => {
    class      => "CQSMappedTable",
    perform    => 1,
    target_dir => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA_table",
    option     => "",
    source_ref => [ "tRNA_1mm_count", ".xml" ],
    groups     => $def->{groups},
    cqs_tools  => $cqstools,
    prefix     => "tRNA_1mm_",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  tRNA_1mm_position => {
    class      => "CQSMappedPosition",
    perform    => 1,
    target_dir => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_tRNA_position",
    option     => "-o " . $def->{task_name} . "_tRNA.position",
    source_ref => "tRNA_1mm_count",
    cqs_tools  => $cqstools,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  smallRNA_1mm_count => {
    class           => "CQSMappedCount",
    perform         => 1,
    target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA",
    option          => $trnacount_option,
    source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
    fastq_files_ref => "identical",
    seqcount_ref    => [ "identical", ".dupcount\$" ],
    cqs_tools       => $cqstools,
    gff_file        => $def->{smallrna_coordinate},
    samtools        => $samtools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  smallRNA_1mm_category => {
    class           => "CQSSmallRNACategory",
    perform         => 1,
    target_dir      => "${cur_target_dir}/topN_bowtie1_genome_cutadapt_1mm_count_smallRNA_category",
    option          => "",
    source_ref      => [ "smallRNA_1mm_count", ".mapped.xml\$" ],
    mirna_count_ref => [ "mirna_1mm_count", ".mapped.xml\$" ],
    groups          => $def->{groups},
    cqs_tools       => $cqstools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  overall => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${cur_target_dir}/overall",
    option     => "",
    source     => {
      individual => [
        "cutadapt", "fastqlen", "identical",
        "bowtie1_genome_cutadapt_topN_1mm_notidentical",
        "bowtie1_genome_cutadapt_topN_1mm",
        "mirna_1mm_count", "miRNA_1mm_count_overlap", "tRNA_1mm_count", "smallRNA_1mm_count",
      ],
      summary => [ "miRNA_1mm_table", "tRNA_1mm_table", "smallRNA_1mm_category", "miRNA_1mm_overlap_position", "tRNA_1mm_position" ],
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
