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
      addPairedFastqToProcessedBam
      addUmiReadsToProcessedBam
      addMutect2
      addSomaticCNV
      addHaplotypecaller
      addCollectAllelicCounts
      addEncodeATACseq
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
    output_file_ext => ".bam",
    output_other_ext => ".bai",
    pbs=> {
      "nodes"     => "1:ppn=8",
      "walltime"  => "2",
      "mem"       => "40gb"
    },
  };
  push @$individual, $task;
  return ($task);
}

sub addPairedFastqToProcessedBam {
  my ($config, $def, $individual, $target_dir, $files_ref) = @_;
  
  my $server_key = getValue($def, "wdl_key", "local");
  my $pipeline_key = "paired_fastq_to_processed_bam";
  my $wdl = $def->{"wdl"};
  my $server = $wdl->{$server_key};
  my $pipeline = $server->{$pipeline_key};

  my $PreProcessing_DoMarkDuplicates = getValue($def, "PreProcessing_DoMarkDuplicates", "false");
  my $genomeForOutputExt=getValue($def, "annovar_buildver", "hg38");

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
      "PreProcessingForVariantDiscovery_GATK4.sample_name" => "SAMPLE_NAME",
      "PreProcessingForVariantDiscovery_GATK4.DoMarkDuplicates" =>$PreProcessing_DoMarkDuplicates,
    },
    "input_list" => {
      #"PreProcessingForVariantDiscovery_GATK4.flowcell_unmapped_bams_list_ref" => [$files_ref,".fastq"]
      "PreProcessingForVariantDiscovery_GATK4.flowcell_unmapped_bams_list_ref" => $files_ref
    },
    output_file_ext => ".".$genomeForOutputExt.".bam",
    output_other_ext => ".".$genomeForOutputExt.".bai",
    pbs=> {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => getValue($def, "PreProcessingForVariantDiscovery_GATK4.memory", "40gb")
    },
  };

  print(Dumper($config->{$task}));
  push @$individual, $task;
  return ($task);
}

sub addHaplotypecaller {
  my ($config, $def, $individual, $target_dir, $files_ref) = @_;
  
  my $server_key = getValue($def, "wdl_key", "local");
  my $pipeline_key = "haplotypecaller";
  my $wdl = $def->{"wdl"};
  my $server = $wdl->{$server_key};
  my $pipeline = $server->{$pipeline_key};

  my $task = $pipeline_key;
  $config->{$task} = {     
    "class" => "CQS::Wdl",
    "target_dir" => "${target_dir}/$pipeline_key",
    "source_ref" => [$files_ref],
    "singularity_image_files_ref" => ["singularity_image_files"],
    "cromwell_jar" => $wdl->{"cromwell_jar"},
    "input_option_file" => $wdl->{"cromwell_option_file"},
    "cromwell_config_file" => $server->{"cromwell_config_file"},
    "wdl_file" => $pipeline->{"wdl_file"},
    "input_json_file" => $pipeline->{"input_file"},
    "input_parameters" => {
      "HaplotypeCallerGvcf_GATK4.input_bam_ref" =>  [$files_ref,".bam\$"],
      "HaplotypeCallerGvcf_GATK4.input_bam_index_ref" =>  [$files_ref,".bai\$"]
    },
    output_file_ext => ".gvcf",
    pbs=> {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => getValue($def, "HaplotypeCallerGvcf_GATK4.memory", "40gb")
    },
  };

  print(Dumper($config->{$task}));
  push @$individual, $task;
  return ($task);
}

sub addUmiReadsToProcessedBam {
  my ($config, $def, $individual, $target_dir, $files_ref, $unmapped_bam_ref) = @_;
  
  my $server_key = getValue($def, "wdl_key", "local");
  my $pipeline_key = "paired_fastq_to_processed_bam";
  my $wdl = $def->{"wdl"};
  my $server = $wdl->{$server_key};
  my $pipeline = $server->{$pipeline_key};
  my $genomeForOutputExt=getValue($def, "annovar_buildver", "hg38");

  my $task = $pipeline_key;
  $config->{$task} = {     
    "class" => "CQS::Wdl",
    "target_dir" => "${target_dir}/$pipeline_key",
    "source_ref" => $files_ref,
    "singularity_image_files_ref" => ["singularity_image_files"],
    "cromwell_jar" => $wdl->{"cromwell_jar"},
    "input_option_file" => $wdl->{"cromwell_option_file"},
    "cromwell_config_file" => $server->{"cromwell_config_file"},
    "wdl_file" => "/home/zhaos/source/perl_cqs/workflow/gatk4-data-processing/processing-for-variant-discovery-gatk4-UMI.wdl",
    "input_json_file" => $pipeline->{"input_file"},
    "input_parameters" => {
      "PreProcessingForVariantDiscovery_GATK4.sample_name" => "SAMPLE_NAME",
      "PreProcessingForVariantDiscovery_GATK4.unmapped_bam_file_ref" => [$unmapped_bam_ref,".bam"],
    },
    "input_list" => {
      "PreProcessingForVariantDiscovery_GATK4.flowcell_unmapped_bams_list_ref" => [$files_ref,".fastq"]
    },
    output_file_ext => ".".$genomeForOutputExt.".bam",
    output_other_ext => ".".$genomeForOutputExt.".bai",
    pbs=> {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => "60gb"
    },
  };
  push @$individual, $task;
  return ($task);
}



sub addMutect2 {
  my ($config, $def, $tasks, $target_dir, $bam_input) = @_;

  my $mutect2_index_dic = {};
  my $mutect2_index_key = "mutect2_Index";
  my $mutect2_prefix = "${bam_input}_muTect2_";

  if (not defined $def->{mutect2_groups}){
    #for wdl mutect2, the result file will use tumor sample name in output
    my $mutect2_groups = {};
    my $groups = $def->{groups};
    if(not defined $groups){
      $groups = $config->{groups};
    }

    if(not defined $groups){
      die "Define groups or mutect2_groups in user definition or configuration";
    }
    
    for my $group_name (sort keys %$groups) {
      my $samples = $groups->{$group_name};
      if (scalar(@$samples) == 1) {
        $mutect2_groups->{$samples->[0]} = $samples;
      }else{
        $mutect2_groups->{$samples->[1]} = $samples;
      }
    }
    $config->{mutect2_groups} = $mutect2_groups;
  }else{
    $config->{mutect2_groups} = $def->{mutect2_groups};
  }
  
  my $mutect2_groups = $config->{mutect2_groups};
  my $mutect2_tumor_only = 0;
  for my $files (values %$mutect2_groups){
    if (scalar(@$files) == 1){
      $mutect2_tumor_only = 1;
      last;
    }
  }

  #die "tumor only " if $mutect2_tumor_only;

  my $mutect2_normal_files=$mutect2_prefix . "_normal_files";
  my $mutect2_tumor_files=$mutect2_prefix . "_tumor_files";
  if($mutect2_tumor_only){
    $config->{$mutect2_tumor_files} = {     
      "class" => "CQS::GroupPickTask",
      "source_ref" => $bam_input,
      "groups_ref" => "mutect2_groups",
      "sample_index_in_group" => 0, 
    };
  }else{
    $config->{$mutect2_normal_files} = {     
      "class" => "CQS::GroupPickTask",
      "source_ref" => $bam_input,
      "groups_ref" => "mutect2_groups",
      "sample_index_in_group" => 0, 
    };
    $config->{$mutect2_tumor_files} = {     
      "class" => "CQS::GroupPickTask",
      "source_ref" => $bam_input,
      "groups_ref" => "mutect2_groups",
      "sample_index_in_group" => 1, 
    };
  }

  my $server_key = getValue($def, "wdl_key", "local");
  my $wdl = $def->{"wdl"};
  my $server = $wdl->{$server_key};
  my $mutect2_pipeline = $server->{"mutect2"};

  my $pon = {};
  if ((not $mutect2_tumor_only) && $mutect2_pipeline->{"perform_mutect2_pon"}){
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
      "source_ref" => [$mutect2_normal_files, ".bam\$"],
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
      output_other_ext => ".vcf.idx",
      pbs=> {
        "nodes"     => "1:ppn=8",
        "walltime"  => "24",
        "mem"       => "70gb"
      },
    };
    push @$tasks, $mutect2_pon;

    $pon = {
      "Mutect2.pon_ref" => [$mutect2_pon, ".vcf\$"],
      "Mutect2.pon_idx_ref" =>[$mutect2_pon, ".vcf.idx\$"],
    };
  }

  my $mutect2_call = $mutect2_prefix . getNextIndex($mutect2_index_dic, $mutect2_index_key) . "_call";
  my $run_funcotator="false";
  if ($def->{ncbi_build} eq "GRCh38") { #based on genome, hg38=true, else false
    $run_funcotator="true";
  }

  my $output_sample_ext;
  if($def->{muTect2_suffix}) {
    $output_sample_ext = $def->{muTect2_suffix};
  }else {
    $output_sample_ext=".hg19";
    if ($def->{ncbi_build} eq "GRCh38") { #based on genome, hg38=true, else false
      $output_sample_ext=".hg38";
    } elsif ($def->{ncbi_build} eq "GRCm38")  {
      $output_sample_ext=".mm10";
    }
  }

  my $output_file_ext;
  my $output_other_ext;
  if ($def->{ncbi_build} eq "GRCh38"){
    $output_file_ext = $output_sample_ext."-filtered.annotated.maf";
    $output_other_ext = $output_sample_ext."-filtered.vcf";
  }else{
    $output_file_ext = $output_sample_ext."-filtered.vcf";
  }

  $config->{$mutect2_call} = {     
    "class" => "CQS::Wdl",
    "target_dir" => "${target_dir}/$mutect2_call",
    "source_ref" => [$mutect2_tumor_files, ".bam\$"],
    "singularity_image_files_ref" => ["singularity_image_files"],
    "cromwell_jar" => $wdl->{"cromwell_jar"},
    "input_option_file" => $wdl->{"cromwell_option_file"},
    "cromwell_config_file" => $server->{"cromwell_config_file"},
    "wdl_file" => $mutect2_pipeline->{"wdl_file"},
    output_file_ext => $output_file_ext,
    output_other_ext => $output_other_ext,
    "input_json_file" => $mutect2_pipeline->{"input_file"},
    "input_parameters" => {
      "Mutect2.intervals" => $def->{covered_bed},
      "Mutect2.ref_fasta" => $def->{ref_fasta},
      "Mutect2.ref_dict" => $def->{ref_fasta_dict},
      "Mutect2.ref_fai" => $def->{ref_fasta} . ".fai",
      "Mutect2.tumor_reads_ref" => [$mutect2_tumor_files, ".bam\$"],
      "Mutect2.tumor_reads_index_ref" => [$mutect2_tumor_files, ".bai\$"],
      "Mutect2.run_funcotator" => $run_funcotator,
    },
    "input_single" => $pon,
    pbs=> {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => "70gb"
    },
  };

  if($mutect2_tumor_only){
    $config->{$mutect2_call}{"input_parameters"}{"Mutect2.normal_reads"} = "";
    $config->{$mutect2_call}{"input_parameters"}{"Mutect2.normal_reads_index"} = "";
  }else{
    $config->{$mutect2_call}{"input_parameters"}{"Mutect2.normal_reads_ref"} = [$mutect2_normal_files, ".bam\$"];
    $config->{$mutect2_call}{"input_parameters"}{"Mutect2.normal_reads_index_ref"} = [$mutect2_normal_files, ".bai\$"];
  }
  push @$tasks, $mutect2_call;
  return ($mutect2_call);
}



sub addSomaticCNV {
  my ($config, $def, $summary, $target_dir, $bam_input) = @_;

  my $somaticCNV_index_dic = {};
  my $somaticCNV_index_key = "somaticCNV_Index";
  my $somaticCNV_prefix = "${bam_input}_somaticCNV_";
  
  my $somaticCNV_normal_files = $somaticCNV_prefix . "_normal_files";
  $config->{$somaticCNV_normal_files} = {     
    "class" => "CQS::GroupPickTask",
    "source_ref" => $bam_input,
    "groups_ref" => "groups",
    "sample_index_in_group" => 0, 
  };
  
  my $somaticCNV_tumor_files = $somaticCNV_prefix . "_tumor_files";
  $config->{$somaticCNV_tumor_files} = {     
    "class" => "CQS::GroupPickTask",
    "source_ref" => $bam_input,
    "groups_ref" => "groups",
    "sample_index_in_group" => 1, 
  };

  my $server_key = getValue($def, "wdl_key", "local");
  my $wdl = $def->{"wdl"};
  my $server = $wdl->{$server_key};
  my $somaticCNV_pipeline = $server->{"somaticCNV"};

  my $pon = {};
  #if ($somaticCNV_pipeline->{"perform_somaticCNV_pon"}){
  if (1) {
    my $pon_pipeline = $server->{"somaticCNV_pon"};

    my $somaticCNV_pon = $somaticCNV_prefix . getNextIndex($somaticCNV_index_dic, $somaticCNV_index_key) . "_pon";
    $config->{$somaticCNV_pon} = {     
      "class" => "CQS::UniqueWdl",
      "target_dir" => "${target_dir}/$somaticCNV_pon",
      "singularity_image_files_ref" => ["singularity_image_files"],
      "cromwell_jar" => $wdl->{"cromwell_jar"},
      "input_option_file" => $wdl->{"cromwell_option_file"},
      "cromwell_config_file" => $server->{"cromwell_config_file"},
      "wdl_file" => $pon_pipeline->{"wdl_file"},
      "input_json_file" => $pon_pipeline->{"input_file"},
      "input_array" => {
        "CNVSomaticPanelWorkflow.normal_bams_ref" => [$somaticCNV_normal_files, ".bam\$"],
        "CNVSomaticPanelWorkflow.normal_bais_ref" => [$somaticCNV_normal_files, ".bai\$"]
      },
      "input_parameters" => {
        "CNVSomaticPanelWorkflow.pon_entity_id" => $config->{general}{task_name},
        "CNVSomaticPanelWorkflow.intervals" => $def->{covered_bed},
      },
      output_file_ext => ".pon.hdf5",
      pbs=> {
        "nodes"     => "1:ppn=8",
        "walltime"  => "48",
        "mem"       => "70gb"
      },
    };
    #push @$summary, $somaticCNV_pon;

#    $pon=[$somaticCNV_pon, ".pon.hdf5\$"];
     $pon = {
       "CNVSomaticPairWorkflow.read_count_pon_ref" => [$somaticCNV_pon, ".pon.hdf5\$"],
# #      "somaticCNV.pon_idx" =>[$somaticCNV_pon, ".vcf.idx\$"],
     };
  }

  my $somaticCNV_call = $somaticCNV_prefix . getNextIndex($somaticCNV_index_dic, $somaticCNV_index_key) . "_call";
  my $run_funcotator="true";
  my $funcotator_ref_version="";
  if ($def->{ncbi_build} eq "GRCh19") { #based on genome, hg38=true, else false
    $funcotator_ref_version="hg19";
  } elsif ($def->{ncbi_build} eq "GRCh38") { #based on genome, hg38=true, else false
    $funcotator_ref_version="hg38";
  } else {
    $run_funcotator="false";
  }
  my $output_genome_ext="hg19";
  if ($def->{ncbi_build} eq "GRCh38") { #based on genome, hg38=true, else false
    $output_genome_ext="hg38";
  } elsif ($def->{ncbi_build} eq "GRCm38")  {
    $output_genome_ext="mm10";
  }

  $config->{$somaticCNV_call} = {     
    "class" => "CQS::Wdl",
    "target_dir" => "${target_dir}/$somaticCNV_call",
    "source_ref" => [$somaticCNV_normal_files, ".bam\$"],
    "singularity_image_files_ref" => ["singularity_image_files"],
    "cromwell_jar" => $wdl->{"cromwell_jar"},
    "input_option_file" => $wdl->{"cromwell_option_file"},
    "cromwell_config_file" => $server->{"cromwell_config_file"},
    "wdl_file" => $somaticCNV_pipeline->{"wdl_file"},
    output_file_ext => ".".$output_genome_ext.".called.seg",
#    output_other_ext => ".".$output_sample_ext."-filtered.vcf",
    "input_json_file" => $somaticCNV_pipeline->{"input_file"},
    "input_parameters" => {
      "CNVSomaticPairWorkflow.intervals" => $def->{covered_bed},
      "CNVSomaticPairWorkflow.funcotator_ref_version" => $funcotator_ref_version,
      "CNVSomaticPairWorkflow.is_run_funcotator" => $run_funcotator,
      "CNVSomaticPairWorkflow.normal_bam_ref" => [$somaticCNV_normal_files, ".bam\$"],
      "CNVSomaticPairWorkflow.normal_bam_idx_ref" => [$somaticCNV_normal_files, ".bai\$"],
      "CNVSomaticPairWorkflow.tumor_bam_ref" => [$somaticCNV_tumor_files, ".bam\$"],
      "CNVSomaticPairWorkflow.tumor_bam_idx_ref" => [$somaticCNV_tumor_files, ".bai\$"],
    },
    "input_single" => $pon,
    pbs=> {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => "70gb"
    },
  };

  if (defined($def->{"CNVSomaticPairWorkflow.common_sites"}) and $def->{"CNVSomaticPairWorkflow.common_sites"} ne "") {
    $config->{$somaticCNV_call}->{input_parameters}->{"CNVSomaticPairWorkflow.common_sites"}=$def->{"CNVSomaticPairWorkflow.common_sites"};
  }

#summary CNV results
  my $somaticCNV_call_summary = $somaticCNV_prefix . getNextIndex($somaticCNV_index_dic, $somaticCNV_index_key) . "_summary";
    $config->{$somaticCNV_call_summary} = {
      class                      => "CQS::UniqueR",
      perform                    => 1,
      target_dir                 => $target_dir . '/' . $somaticCNV_call_summary,
      rtemplate                  => "../CNV/GATKsomaticCNVSummary.R",
      parameterSampleFile1_ref   => [ $somaticCNV_call, ".called.seg\$" ],
 #     parameterFile1_ref         => [ $cnvAnnotationGenesPlot, ".position.txt.slim" ],
#      parameterSampleFile2       => $def->{onco_options},
#      parameterSampleFile3       => $def->{onco_sample_groups},
      output_to_result_directory => 1,
      output_file                => "",
      output_file_ext            => ".allCNV.seg;.allCNV.filter.seg",
      sh_direct                  => 1,
      'pbs'                      => {
        'nodes'    => '1:ppn=1',
        'mem'      => '20gb',
        'walltime' => '10'
      },
    };

  #push @$summary, $somaticCNV_call;
  return ($somaticCNV_call);
}

sub addCollectAllelicCounts {
  my ($config, $def, $individual, $target_dir,$bam_input,$common_sites) = @_;
  
  my $server_key = getValue($def, "wdl_key", "local");
  my $pipeline_key = "CollectAllelicCounts";
  my $wdl = $def->{"wdl"};
  my $server = $wdl->{$server_key};
  my $pipeline = $server->{$pipeline_key};

  my $genomeForOutputExt=getValue($def, "annovar_buildver", "hg38");
  my $task = "${bam_input}_CollectAllelicCounts";

  #only need tumor bam.
  #$config->{$mutect2_tumor_files} was predefined in mutect2 task, here use it to get tumor samples/bams
  #so don't need redefined it again
  #my $mutect2_tumor_files = $mutect_prefix . "_tumor_files";

  # if (not defined $def->{mutect2_groups}){
  #   #for wdl mutect2, the result file will use tumor sample name in output
  #   my $mutect2_groups = {};
  #   my $groups = $def->{groups};
  #   for my $group_name (sort keys %$groups) {
  #     my $samples = $groups->{$group_name};
  #     if (scalar(@$samples) == 1) {
  #       $mutect2_groups->{$samples->[0]} = $samples;
  #     }else{
  #       $mutect2_groups->{$samples->[1]} = $samples;
  #     }
  #   }
  #   $config->{mutect2_groups} = $mutect2_groups;
  # }else{
  #   $config->{mutect2_groups} = $def->{mutect2_groups};
  # }
  # $config->{$mutect2_tumor_files} = {     
  #   "class" => "CQS::GroupPickTask",
  #   "source_ref" => $bam_input,
  #   "groups_ref" => "mutect2_groups",
  #   "sample_index_in_group" => 1, 
  # };

  #if addMutect2 was run before, $mutect2_tumor_files exists in $config. Can use it to skip normal bam files
  my $mutect2_prefix = "${bam_input}_muTect2_";
  my $mutect2_tumor_files = $mutect2_prefix . "_tumor_files";
  if (exists($config->{$mutect2_tumor_files})) {
    $bam_input=$mutect2_tumor_files;
  }

  $config->{$task} = {     
    "class" => "CQS::Wdl",
    "target_dir" => "${target_dir}/$task",
    "source_ref" => $bam_input,
    "singularity_image_files_ref" => ["singularity_image_files"],
    "cromwell_jar" => $wdl->{"cromwell_jar"},
    "input_option_file" => $wdl->{"cromwell_option_file"},
    "cromwell_config_file" => $server->{"cromwell_config_file"},
    "wdl_file" => $pipeline->{"wdl_file"},
    "input_json_file" => $pipeline->{"input_file"},
    "input_parameters" => {
      "CollectAllelicCountsWorkflow.ref_fasta" => $def->{ref_fasta},
      "CollectAllelicCountsWorkflow.ref_fasta_dict" => $def->{ref_fasta_dict},
      "CollectAllelicCountsWorkflow.ref_fasta_fai" => $def->{ref_fasta} . ".fai",
      "CollectAllelicCountsWorkflow.tumor_bam_ref" =>  [$bam_input, ".bam\$"],
      "CollectAllelicCountsWorkflow.tumor_bam_idx_ref" =>  [$bam_input, ".bai\$"],
    },
    "input_single" => {
      "CollectAllelicCountsWorkflow.common_sites_ref" =>  [$common_sites, ".intervals\$"],
    },
    output_file_ext => ".$genomeForOutputExt.allelicCounts.tsv",
#    output_other_ext => ".".$genomeForOutputExt.".bai",
    pbs=> {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => getValue($def, "CollectAllelicCountsWorkflow.memory", "40gb")
    },
  };

  print(Dumper($config->{$task}));
  push @$individual, $task;
  return ($task);
}

sub addEncodeATACseq {
  my ($config, $def, $individual, $target_dir, $files_ref) = @_;

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
  my $pipeline_key = "encode_atacseq";
  my $wdl = $def->{"wdl"};
  my $server = $wdl->{$server_key};
  my $pipeline = $server->{$pipeline_key};

  my $task = $pipeline_key;
  $config->{$task} = {     
    "class" => "CQS::Wdl",
    "option" => "--no-build-singularity",
    "target_dir" => "${target_dir}/$pipeline_key",
    "source_ref" => $files_ref,
    "singularity_image_files_ref" => ["singularity_image_files"],
    "cromwell_jar" => $wdl->{"cromwell_jar"},
    "input_option_file" => $wdl->{"cromwell_option_file"},
    "cromwell_config_file" => $server->{"cromwell_config_file"},
    "wdl_file" => $pipeline->{"wdl_file"},
    "input_json_file" => $pipeline->{"input_file"},
    "input_parameters" => {
      "atac.title" => "SAMPLE_NAME",
      "atac.description" => "SAMPLE_NAME",
      "atac.genome_tsv" => getValue($def, "encode_atacseq_genome_tsv"),
      "atac.paired_end" => is_paired_end($def) ? "true" : "false",
      "atac.adapter" => getValue($def, "adapter", ""),
      "atac.fastqs_rep1_R1_ref" => [$fastq_1],
      "atac.fastqs_rep1_R2_ref" => [$fastq_2]
    },
    output_to_same_folder => 0,
    cromwell_finalOutputs => 0,
    check_output_file_pattern => "metadata.json",
    output_file_ext => "atac/",
    use_caper => 1,
    pbs=> {
      "nodes"     => "1:ppn=8",
      "walltime"  => "12",
      "mem"       => "40gb"
    },
  };
  push @$individual, $task;

  my $croo_task = $task . "_croo";
  $config->{$croo_task} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => "${target_dir}/$croo_task",
    interpretor => "python3",
    program => "../Chipseq/croo.py",
    option => "-n __NAME__ --croo " . getValue($def, "croo", "croo") . " --out_def_json " . getValue($def, "croo_out_def_json"),
    source_arg => "-i",
    source_ref => [$task],
    output_arg => "-o",
    output_file_prefix => "",
    output_file_ext => "__NAME__/peak/overlap_reproducibility/overlap.optimal_peak.narrowPeak.gz",
    output_to_same_folder => 1,
    can_result_be_empty_file => 0,
    sh_direct   => 1,
    pbs => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push @$individual, $croo_task;

  return ($task);
}

1;
