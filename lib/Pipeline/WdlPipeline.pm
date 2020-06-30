#!/usr/bin/perl
package Pipeline::WdlPipeline;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(addMutect2)
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub addMutect2 {
  my ($config, $def, $summary, $target_dir, $bam_input) = @_;

  my $mutect2_index_dic = {};
  my $mutect2_index_key = "mutect2_Index";
  my $mutect2_prefix = "${bam_input}_muTect2_";
  my $wdl_key = getValue($def, "wdl_key", "local");
  
  my $mutect2_normal_files = $mutect2_prefix . "_normal_files";
  $config->{$mutect2_normal_files} = {     
    "class" => "CQS::GroupPickTask",
    "source_ref" => [$bam_input],
    "groups_ref" => ["groups"],
    "sample_index_in_group" => 0, 
  };
  
  my $mutect2_tumor_files = $mutect2_prefix . "_tumor_files";
  $config->{$mutect2_tumor_files} = {     
    "class" => "CQS::GroupPickTask",
    "source_ref" => [$bam_input],
    "groups_ref" => ["groups"],
    "sample_index_in_group" => 1, 
  };

  my $mutect2_def = $def->{wdl}{$wdl_key}{mutect2};

  my $pon = {};
  if ($mutect2_def->{perform_mutect2_pon}){
    my $pon_def = $def->{wdl}{$wdl_key}{mutect2_pon};
    my $mutect2_pon = $mutect2_prefix . getNextIndex($mutect2_index_dic, $mutect2_index_key) . "_pon";
    $config->{$mutect2_pon} = {     
      "class" => "CQS::UniqueWdl",
      "target_dir" => "${target_dir}/$mutect2_pon",
      "cromwell_config_file" => $def->{cromwell_config_file}{$wdl_key},
      "cromwell_jar" => getValue($def, "cromwell_jar"),
      "input_option_file" => getValue($def, "cromwell_option_file"),
      "singularity_image_files_ref" => ["singularity_image_files"],
      "wdl_file" => $pon_def->{wdl_file},
      "input_json_file" => $pon_def->{input_file},
      "input_array" => {
        "Mutect2_Panel.normal_bams_ref" => [$mutect2_normal_files, ".bam\$"],
        "Mutect2_Panel.normal_bais_ref" => [$mutect2_normal_files, ".bai\$"]
      },
      "input_parameters" => {
        "Mutect2_Panel.pon_name" => $config->{general}{task_name},
        "Mutect2_Panel.intervals" => $def->{covered_bed}
      },
      output_file_ext => ".vcf",
      output_other_exts => ".vcf.idx",
      pbs=> {
        "nodes"     => "1:ppn=8",
        "walltime"  => "2",
        "mem"       => "40gb"
      },
    };
    push @$summary, $mutect2_pon;

    $pon = {
      "Mutect2.pon" => [$mutect2_pon, ".vcf\$"],
      "Mutect2.pon_idx" =>[$mutect2_pon, ".vcf.idx\$"],
    };
  }

  my $mutect2_call = $mutect2_prefix . getNextIndex($mutect2_index_dic, $mutect2_index_key) . "_call";
  $config->{$mutect2_call} = {     
    "class" => "CQS::Wdl",
    "target_dir" => "${target_dir}/$mutect2_call",
    "source_ref" => [$mutect2_normal_files, ".bam\$"],
    "cromwell_config_file" => $def->{cromwell_config_file}{$wdl_key},
    "cromwell_jar" => getValue($def, "cromwell_jar"),
    "singularity_image_files_ref" => ["singularity_image_files"],
    "wdl_file" => $def->{wdl}{$wdl_key}{mutect2}{wdl_file},
    "input_json_file" => $def->{wdl}{$wdl_key}{mutect2}{input_file},
    "input_option_file" => getValue($def, "cromwell_option_file"),
    "input_parameters" => {
      "Mutect2.intervals" => $def->{covered_bed},
      "Mutect2.ref_fasta" => $def->{ref_fasta},
      "Mutect2.ref_dict" => $def->{ref_fasta_dict},
      "Mutect2.ref_fai" => $def->{ref_fasta} . ".fai",
      "Mutect2.normal_reads_ref" => [$mutect2_normal_files, ".bam\$"],
      "Mutect2.normal_reads_index_ref" => [$mutect2_normal_files, ".bai\$"],
      "Mutect2.tumor_reads_ref" => [$mutect2_tumor_files, ".bam\$"],
      "Mutect2.tumor_reads_index_ref" => [$mutect2_tumor_files, ".bai\$"],
    },
    "input_single" => $pon,
    pbs=> {
      "nodes"     => "1:ppn=8",
      "walltime"  => "2",
      "mem"       => "40gb"
    },
  };
  push @$summary, $mutect2_call;
  return ($mutect2_call);
}