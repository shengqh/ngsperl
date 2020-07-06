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
use POSIX qw(strftime);

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(
      addPairedFastqToUnmappedBam
      addMutect2
    )
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub addPairedFastqToUnmappedBam {
  my ($config, $def, $individual, $target_dir, $files_ref) = @_;

  my $datestring = strftime("%Y-%m-%dT%H:%M:%S%z", localtime);

  my $fastq_1 = "fastq_1";
  $config->{$fastq_1} = {     
    "class" => "CQS::FilePickTask",
    "source_ref" => $files_ref,
    "sample_index" => 0, 
  };
  
  my $fastq_2 = "fastq_2";
  $config->{$fastq_2} = {     
    "class" => "CQS::FilePickTask",
    "source_ref" => $files_ref,
    "sample_index" => 1, 
  };
  
  my $server_key = getValue($def, "wdl_key", "local");
  my $pipeline_key = "paired_fastq_to_unmapped_bam";
  my $wdl = $def->{"wdl"};
  my $server = $wdl->{$server_key};
  my $pipeline = $server->{$pipeline_key};

  my $task = $pipeline_key;
  $config->{$task} = {     
    "class" => "CQS::Wdl",
    "target_dir" => "${target_dir}/$pipeline_key",
    "source_ref" => $files_ref,
    "singularity_image_files_ref" => ["singularity_image_files"],
    "cromwell_jar" => $wdl->{"cromwell_jar"},
    "input_option_file" => $wdl->{"cromwell_option_file"},
    "cromwell_config_file" => $server->{"cromwell_config_file"},
    "wdl_file" => $pipeline->{"wdl_file"},
    "input_json_file" => $pipeline->{"input_file"},
    "input_parameters" => {
      "ConvertPairedFastQsToUnmappedBamWf.readgroup_name" => "SAMPLE_NAME",
      "ConvertPairedFastQsToUnmappedBamWf.sample_name" => "SAMPLE_NAME",
      "ConvertPairedFastQsToUnmappedBamWf.library_name" => "SAMPLE_NAME",
      "ConvertPairedFastQsToUnmappedBamWf.platform_unit" => "illumina",
      "ConvertPairedFastQsToUnmappedBamWf.run_date" => $datestring,
      "ConvertPairedFastQsToUnmappedBamWf.platform_name" => "illumina",
      "ConvertPairedFastQsToUnmappedBamWf.sequencing_center" => "Unknown",
      "ConvertPairedFastQsToUnmappedBamWf.fastq_1_ref" => [$fastq_1],
      "ConvertPairedFastQsToUnmappedBamWf.fastq_2_ref" => [$fastq_2]
    },
    pbs=> {
      "nodes"     => "1:ppn=8",
      "walltime"  => "2",
      "mem"       => "40gb"
    },
  };
  push @$individual, $task;
  return ($task);
}

sub addMutect2 {
  my ($config, $def, $summary, $target_dir, $bam_input) = @_;

  my $mutect2_index_dic = {};
  my $mutect2_index_key = "mutect2_Index";
  my $mutect2_prefix = "${bam_input}_muTect2_";
  
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

  my $server_key = getValue($def, "wdl_key", "local");
  my $wdl = $def->{"wdl"};
  my $server = $wdl->{$server_key};
  my $mutect2_pipeline = $server->{"mutect2"};

  my $pon = {};
  if ($mutect2_pipeline->{"perform_mutect2_pon"}){
    my $pon_pipeline = $server->{"mutect2_pon"};

    my $mutect2_pon = $mutect2_prefix . getNextIndex($mutect2_index_dic, $mutect2_index_key) . "_pon";
    $config->{$mutect2_pon} = {     
      "class" => "CQS::UniqueWdl",
      "target_dir" => "${target_dir}/$mutect2_pon",
      "singularity_image_files_ref" => ["singularity_image_files"],
      "cromwell_jar" => $wdl->{"cromwell_jar"},
      "input_option_file" => $wdl->{"cromwell_option_file"},
      "cromwell_config_file" => $server->{"cromwell_config_file"},
      "wdl_file" => $pon_pipeline->{"wdl_file"},
      "input_json_file" => $pon_pipeline->{"input_file"},
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
    "singularity_image_files_ref" => ["singularity_image_files"],
    "cromwell_jar" => $wdl->{"cromwell_jar"},
    "input_option_file" => $wdl->{"cromwell_option_file"},
    "cromwell_config_file" => $server->{"cromwell_config_file"},
    "wdl_file" => $mutect2_pipeline->{"wdl_file"},
    "input_json_file" => $mutect2_pipeline->{"input_file"},
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

1;