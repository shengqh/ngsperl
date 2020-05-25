#!/usr/bin/perl
package CQS::TestWdl;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::ClassFactory;
use CQS::Wdl;

my $config = {
  general => {
    task_name => "wdl",
    email => "quanhu.sheng.1\@vumc.org",
  },
  normal_files => {
    "S1" => ["NormalFile1.bam", "NormalFile1.bam.bai"],
    "S2" => ["NormalFile2.bam", "NormalFile2.bam.bai"],
  },
  tumor_files => {
    "S1" => ["TumorFile1.bam", "TumorFile1.bam.bai"],
    "S2" => ["TumorFile2.bam", "TumorFile2.bam.bai"],
  },
  "test" => {
    "class" => "CQS::Wdl",
    "target_dir" => "c:/temp",
    "source_ref" => ["normal_files", ".bam\$"],
    "cromwell_config_file" => "C:/Users/sheng/git/ngsperl/data/wdl.input.json",
    "cromwell_jar" => "C:/Users/sheng/git/ngsperl/data/wdl.input.json",
    "pipeline_file" => "C:/Users/sheng/git/ngsperl/data/wdl.input.json",
    "wdl_file" => "C:/Users/sheng/git/ngsperl/data/wdl.input.json",
    "input_option_file" => "C:/Users/sheng/git/ngsperl/data/wdl.input.json",
    "input_json_file" => "C:/Users/sheng/git/ngsperl/data/wdl.input.json",
    "input_parameters" => {
      "Mutect2.normal_reads_ref" => ["normal_files", ".bam\$"],
      "Mutect2.normal_reads_index_ref" => ["normal_files", ".bam.bai"],
      "Mutect2.tumor_reads_ref" => ["tumor_files", ".bam\$"],
      "Mutect2.tumor_reads_index_ref" => ["tumor_files", ".bam.bai"],
      "Mutect2.pon" => "panel_of_normal.vcf",
      "Mutect2.pon_idx" =>"panel_of_normal.vcf.idx"
    },
    pbs                       => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "12",
      "mem"       => "40gb"
    },
  },
};

performTask($config, "test");

1