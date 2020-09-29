#!/usr/bin/perl
package CQS::TestDESeq2config;

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

my $config = {
  general => {
    task_name => "wdl",
    email => "quanhu.sheng.1\@vumc.org",
    covariance_file => "E:/sqh/programs/perl/ngsperl/data/covariance.txt", 
  },
  groups => {
    Normal => ["Normal"],
    Tumor => ["Tumor"],
  },
  pairs => {
    "Tumor_vs_Normal" => {
      groups => ["Normal", "Tumor"],
      covariances => ["Age", "Gender"],
      formula => "~Age+Gender+Condition"
    }
  },
  "deseq2" => {
    "class" => "Comparison::DESeq2covariance",
    "target_dir" => "E:/temp",
    "source_ref" => "pairs",
    "groups_ref" => "groups",
    "countfile" => "E:/sqh/programs/perl/ngsperl/data/covariance.count.txt",
    pbs                       => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "12",
      "mem"       => "40gb"
    },
  },
};

performTask($config, "deseq2");

1