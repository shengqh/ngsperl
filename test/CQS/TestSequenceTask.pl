#!/usr/bin/env perl
use strict;
use warnings;

use Test::More tests => 2;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use CQS::SequenceTaskSlurmSlim;

my $annovar_protocol  = "refGene,snp138,cosmic70";
my $annovar_operation = "g,f,f";
my $annovar_param     = "-protocol ${annovar_protocol} -operation ${annovar_operation} --remove";
my $annovar_db        = "/scratch/cqs/shengq2/references/annovar/humandb/";
my $email             = "quanhu.sheng.1\@vanderbilt.edu";

my $config = {
  general => {
    task_name => "annovar",
    cluster   => "slurm"
  },
  annovar => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "e:/temp/annovar",
    option     => $annovar_param,
    source     => {
      "Binarypheno1" => ["/scratch/cqs/shengq2/evan/20170823_evan_annovar/vcf/ExomeChipBinarypheno1.vcf"],
      "RVSPLinear"   => ["/scratch/cqs/shengq2/evan/20170823_evan_annovar/vcf/ExomeChipRVSPLinear.vcf"],
    },
    annovar_db => $annovar_db,
    buildver   => "hg19",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  GTEx => {
    class       => "eQTL::GTEx",
    perform     => 1,
    target_dir  => "e:/temp/GTEx",
    option      => "",
    source_ref  => [ "annovar", ".tsv\$" ],
    GTEx_folder => "/scratch/cqs/shengq2/references/GTEx",
    sh_direct   => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  WebGestaltR => {
    class        => "Annotation::WebGestaltR",
    perform      => 1,
    target_dir   => "e:/temp/WebGestaltR",
    option       => "",
    source_ref   => [ "GTEx", ".tsv.genes\$" ],
    annoavar_ref => [ "annovar", ".tsv\$" ],
    organism     => "hsapiens",
    sh_direct    => 1,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTaskSlurm2",
    perform    => 1,
    target_dir => "e:/temp/sequencetask",
    option     => "",
    source     => {
      "step1" => [ "annovar", "GTEx", "WebGestaltR" ]
    },
    sh_direct => 1,
    cluster   => "slurm",
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  }
};

my $refpbsmap = get_ref_section_pbs( $config, "WebGestaltR", "source" );
is_deeply(
  $refpbsmap,
  {
    Binarypheno1 => ['e:/temp/GTEx/pbs/Binarypheno1_gt.pbs'],
    RVSPLinear   => ['e:/temp/GTEx/pbs/RVSPLinear_gt.pbs']
  }
);

my $myclass = instantiate("CQS::SequenceTask");
my $deppbsmap = $myclass->get_dependent_pbs_map( $config, "sequencetask" );
is_deeply(
  $deppbsmap,
  {
    annovar => {},
    GTEx    => {
      Binarypheno1 => {
        'e:/temp/annovar/pbs/Binarypheno1_ann.pbs' => 1,
      },
      RVSPLinear => {
        'e:/temp/annovar/pbs/RVSPLinear_ann.pbs' => 1,
      }
    },
    WebGestaltR => {
      Binarypheno1 => {
        'e:/temp/annovar/pbs/Binarypheno1_ann.pbs' => 1,
        'e:/temp/GTEx/pbs/Binarypheno1_gt.pbs'     => 1
      },
      RVSPLinear => {
        'e:/temp/annovar/pbs/RVSPLinear_ann.pbs' => 1,
        'e:/temp/GTEx/pbs/RVSPLinear_gt.pbs'     => 1
      },
    },
  }
);

performTask($config, 'sequencetask');

1;
