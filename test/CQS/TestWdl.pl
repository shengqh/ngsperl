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
    "HpEtHOBA1_Ear" => ["/scratch/cqs/zhaos/WilsonKeith/20200330_EtHOBAMouseSomatic/20200330_preprocessing_cromwell/result/cromwell_finalOutputs/HpEtHOBA1_Ear.mm10.bam", 
                        "/scratch/cqs/zhaos/WilsonKeith/20200330_EtHOBAMouseSomatic/20200330_preprocessing_cromwell/result/cromwell_finalOutputs/HpEtHOBA1_Ear.mm10.bam.bai"],
  },
  tumor_files => {
    "HpEtHOBA1_Ear" => ["/scratch/cqs/zhaos/WilsonKeith/20200330_EtHOBAMouseSomatic/20200330_preprocessing_cromwell/result/cromwell_finalOutputs/HpEtHOBA1_Stomach.mm10.bam",
                        "/scratch/cqs/zhaos/WilsonKeith/20200330_EtHOBAMouseSomatic/20200330_preprocessing_cromwell/result/cromwell_finalOutputs/HpEtHOBA1_Stomach.mm10.bai"],
  },
  "test" => {
    "class" => "CQS::Wdl",
    "target_dir" => "/scratch/cqs/shengq2/temp",
    "source_ref" => ["normal_files", ".bam\$"],
    "cromwell_config_file" => "/home/zhaos/source/perl_cqs/test/cromwell/cromwell.examples.local.conf",
    "cromwell_jar" => "/scratch/cqs/zhaos/test/cromwell/cromwell-47.jar",
    "wdl_file" => "/home/zhaos/source/perl_cqs/test/cromwell/processing-for-variant-discovery-gatk4-fromPairEndFastq.wdl",
    "input_option_file" => "/home/zhaos/source/perl_cqs/workflow/cromwell.options.json",
    "input_json_file" => "/home/shengq2/program/ngsperl/data/wdl.input.json",
    "input_parameters" => {
      "Mutect2.normal_reads_ref" => ["normal_files", ".bam\$"],
      "Mutect2.normal_reads_index_ref" => ["normal_files", ".bai\$"],
      "Mutect2.tumor_reads_ref" => ["tumor_files", ".bam\$"],
      "Mutect2.tumor_reads_index_ref" => ["tumor_files", ".bai\$"],
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