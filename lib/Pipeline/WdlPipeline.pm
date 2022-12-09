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
use List::Util qw(max);
use POSIX qw(strftime);

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(
      addPairedFastqToUnmappedBam
      addPairedFastqToProcessedBam
      addUmiReadsToProcessedBam
      is_muTect2_tumor_only
      addMutect2Wdl
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
      "PreProcessingForVariantDiscovery_GATK4.SamToFastqAndBwaMem.num_cpu" => "8",
      "PreProcessingForVariantDiscovery_GATK4.DoMarkDuplicates" => $PreProcessing_DoMarkDuplicates,
    },
    "input_list" => {
      #"PreProcessingForVariantDiscovery_GATK4.flowcell_unmapped_bams_list_ref" => [$files_ref,".fastq"]
      "PreProcessingForVariantDiscovery_GATK4.flowcell_unmapped_bams_list_ref" => $files_ref
    },
    output_file_ext => ".".$genomeForOutputExt.".bam",
    output_other_ext => ".".$genomeForOutputExt.".bai",
    pbs=> {
      "nodes"     => "1:ppn=8",
      "walltime"  => "72",
      "mem"       => getValue($def, "PreProcessingForVariantDiscovery_GATK4.memory", "70gb")
    },
  };

  #print(Dumper($config->{$task}));
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

sub is_muTect2_tumor_only {
  my $config = shift;    
  my $mutect2_groups = $config->{mutect2_groups};
  #print(Dumper($mutect2_groups));
  for my $files (values %$mutect2_groups){
    if (scalar(@$files) == 1){
      return(1);
    }
  }
  return(0);
}

sub addMutect2Wdl {
  my ($config, $def, $tasks, $target_dir, $mutect2_prefix, $mutect2_option, $is_pon, $mutect2_tumor_only, $mutect2_tumor_files, $mutect2_normal_files, $pon, $pon_idx) = @_;

  # if(!$is_pon){
  #   print("tumors=" . $mutect2_tumor_files . "\n");
  #   if (!$mutect2_tumor_only){
  #     print("normal=" . $mutect2_normal_files . "\n");
  #   }
  # }

  my $server_key = getValue($def, "wdl_key", "local");
  my $wdl = $def->{"wdl"};
  my $server = $wdl->{$server_key};
  my $mutect2_pipeline = $server->{"mutect2"};

  my $mutect2_call = $mutect2_prefix . getNextIndex($def, $mutect2_prefix) . "_call_wdl";
  my $run_funcotator="false";
  if ($def->{ncbi_build} eq "GRCh38") { #based on genome, hg38=true, else false
    $run_funcotator="true";
  }

  my $run_orientation_bias_mixture_model_filter = getValue($def, "Mutect2.run_orientation_bias_mixture_model_filter", "true");

  my $output_sample_ext = "";
  # if($def->{muTect2_suffix}) {
  #   $output_sample_ext = $def->{muTect2_suffix};
  # }else {
  #   $output_sample_ext=".hg19";
  #   if ($def->{ncbi_build} eq "GRCh38") { #based on genome, hg38=true, else false
  #     $output_sample_ext=".hg38";
  #   } elsif ($def->{ncbi_build} eq "GRCm38")  {
  #     $output_sample_ext=".mm10";
  #   }
  # }

  my $output_file_ext;
  my $output_other_ext;
  if ( $def->{ncbi_build} eq "GRCh38" && ! $is_pon ){
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
    "output_file_ext" => $output_file_ext,
    "output_other_ext" => $output_other_ext,
    "input_json_file" => $mutect2_pipeline->{"input_file"},
    "use_filename_in_result" => 1,
    "input_parameters" => {
      "Mutect2.m2_extra_args" => $mutect2_option,
      "Mutect2.intervals" => $def->{covered_bed},
      "Mutect2.ref_fasta" => $def->{ref_fasta},
      "Mutect2.ref_dict" => $def->{ref_fasta_dict},
      "Mutect2.ref_fai" => $def->{ref_fasta} . ".fai",
      "Mutect2.tumor_reads_ref" => [$mutect2_tumor_files, ".bam\$"],
      "Mutect2.tumor_reads_index_ref" => [$mutect2_tumor_files, ".bai\$"],
      "Mutect2.run_funcotator" => $run_funcotator,
      "Mutect2.run_orientation_bias_mixture_model_filter" => $run_orientation_bias_mixture_model_filter
    },
    "input_single" => {},
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

  if($is_pon){
    $config->{$mutect2_call}{"input_parameters"}{"Mutect2.pon"} = "";
    $config->{$mutect2_call}{"input_parameters"}{"Mutect2.pon_idx"} = "";
    $config->{$mutect2_call}{"input_parameters"}{"Mutect2.gnomad"} = "";
    $config->{$mutect2_call}{"input_parameters"}{"Mutect2.gnomad_idx"} = "";
    $config->{$mutect2_call}{"input_parameters"}{"Mutect2.variants_for_contamination"} = "";
    $config->{$mutect2_call}{"input_parameters"}{"Mutect2.variants_for_contamination_idx"} = "";
    $config->{$mutect2_call}{"input_parameters"}{"Mutect2.run_funcotator"} = "false";
    $config->{$mutect2_call}{"input_parameters"}{"Mutect2.funco_reference_version"} = "";
    $config->{$mutect2_call}{"input_parameters"}{"Mutect2.funco_data_sources_tar_gz"} = "";
    $config->{$mutect2_call}{"input_parameters"}{"Mutect2.funco_transcript_selection_list"} = "";
  }else{
    if( -e $pon ){
      $config->{$mutect2_call}{"input_parameters"}{"Mutect2.pon"} = $pon;
      $config->{$mutect2_call}{"input_parameters"}{"Mutect2.pon_idx"} = $pon_idx;
    }else{
      if(defined $pon){
        if ((index($pon, '/') != -1) || (index($pon, '/') != -1)){
          die "file not exists: " . $pon
        }
        if ((index($pon_idx, '/') != -1) || (index($pon_idx, '/') != -1)){
          die "file not exists: " . $pon_idx
        }
        $config->{$mutect2_call}{"input_single"}{"Mutect2.pon_ref"} = $pon;
        $config->{$mutect2_call}{"input_single"}{"Mutect2.pon_idx_ref"} = $pon_idx;
      }else{
        $config->{$mutect2_call}{"input_parameters"}{"Mutect2.pon"} = "";
        $config->{$mutect2_call}{"input_parameters"}{"Mutect2.pon_idx"} = "";
      }
    }
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
  if ($def->{ncbi_build} eq "GRCh19" or $def->{ncbi_build} eq "GRCh37") { #based on genome, hg38=true, else false
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

  my $somaticCNV_call_output_ext=".called.seg";
  if ($run_funcotator eq "true") {
    $somaticCNV_call_output_ext=$somaticCNV_call_output_ext."; .called.seg.funcotated.tsv";
  }
  $config->{$somaticCNV_call} = {     
    "class" => "CQS::Wdl",
    "target_dir" => "${target_dir}/$somaticCNV_call",
    "source_ref" => [$somaticCNV_tumor_files, ".bam\$"],
    "singularity_image_files_ref" => ["singularity_image_files"],
    "cromwell_jar" => $wdl->{"cromwell_jar"},
    "input_option_file" => $wdl->{"cromwell_option_file"},
    "cromwell_config_file" => $server->{"cromwell_config_file"},
    "wdl_file" => $somaticCNV_pipeline->{"wdl_file"},
    "use_filename_in_result" => 1,
    output_file_ext => $somaticCNV_call_output_ext,
#    output_file_ext => ".".$output_genome_ext.".called.seg",
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
      parameterSampleFile2_ref   => [ $somaticCNV_call, ".called.seg.funcotated.tsv\$" ],
#      parameterFile1_ref         => [ $cnvAnnotationGenesPlot, ".position.txt.slim" ],
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
  my ($config, $def, $individual, $target_dir, $files_ref, $task) = @_;

  my $server_key = getValue($def, "wdl_key", "local");
  my $pipeline_key = "encode_atacseq";
  my $wdl = $def->{"wdl"};
  my $server = $wdl->{$server_key};
  my $pipeline = $server->{$pipeline_key};

  if(not defined $task){
    $task = $pipeline_key;
  }

  #for encode atac adapter is required.
  my $adapter = getValue($def, "adapter");
  #my $adapter = getValue($def, "perform_cutadapt", 0) ? getValue($def, "adapter", "") : "";
  #print("adapter = " . $adapter . "\n");
  my $is_paired_end = is_paired_end($def);
  my $encode_option = getValue($def, "encode_option", "");
  my $folder_suffix = $encode_option =~ /slurm/ ? "_slurm" : "";
  my $sh_direct = $encode_option =~ /slurm/ ? 1 : 0;

  my $task_folder = "${task}${folder_suffix}";

  my $encode_atac_inputs = getValue($def, "encode_atac_inputs", {});

  my $encode_cpu = getValue($def, "encode_atac_cpu", "16");

  $config->{$task} = {     
    "class" => "CQS::Wdl",
    #"option" => "--no-build-singularity",
    "option" => $encode_option,
    "target_dir" => "${target_dir}/${task_folder}",
    "singularity_image_files_ref" => ["singularity_image_files"],
    "cromwell_jar" => $wdl->{"cromwell_jar"},
    "input_option_file" => $wdl->{"cromwell_option_file"},
    "cromwell_config_file" => $server->{"cromwell_config_file"},
    "wdl_file" => $pipeline->{"wdl_file"},
    "input_json_file" => $pipeline->{"input_file"},
    "input_parameters" => merge_hash_left_precedent($encode_atac_inputs, {
      "atac.title" => "SAMPLE_NAME",
      "atac.description" => "SAMPLE_NAME",
      "atac.genome_tsv" => getValue($def, "encode_atacseq_genome_tsv"),
      "atac.paired_end" => $is_paired_end ? "true" : "false",
      "atac.adapter" => $adapter,
      "atac.align_cpu" => $encode_cpu,
      "atac.align_mem_factor" => 2,
      "atac.filter_cpu" => $encode_cpu,
      "atac.filter_mem_factor" => 2,
      "atac.bam2ta_mem_factor" => 1,
    }),
    output_to_same_folder => 0,
    cromwell_finalOutputs => 0,
    check_output_file_pattern => "metadata.json",
    output_file_ext => "atac/",
    use_caper => 1,
    sh_direct   => $sh_direct,
    pbs=> {
      "nodes"     => "1:ppn=$encode_cpu",
      "walltime"  => getValue($def, "encode_atac_walltime", "24"),
      "mem"       => getValue($def, "encode_atac_men", "40gb"),
    },
  };

  my $input_parameters_is_vector = {};
  if(defined $def->{replicates}){
    my $replicates = $def->{replicates};
    $config->{$task}{"source"} = $replicates;
    $config->{$task}{"input_parameters"}{"atac.true_rep_only"} = getValue($def, "atac.true_rep_only", "true");

    my $max_len = 0;
    for my $values (values %$replicates){
      $max_len = max($max_len, scalar(@$values))
    }
    my $pick_index = 0;
    while($pick_index < $max_len){
      my $pick_str = $pick_index + 1;
      my $group_i = "group_" . $pick_str;
      $config->{$group_i} = {     
        "class" => "CQS::GroupPickTask",
        "source_ref" => $files_ref,
        "groups"     => $def->{replicates},
        "sample_index_in_group" => $pick_index,
      };

      my $fastq_1 = "fastq_${pick_str}_1";
      $config->{$fastq_1} = {     
        "class" => "CQS::FilePickTask",
        "source_ref" => [$group_i],
        "sample_index" => 0, 
      };
      $config->{$task}{input_parameters}{"atac.fastqs_rep${pick_str}_R1_ref"} = [$fastq_1];
      $input_parameters_is_vector->{"atac.fastqs_rep${pick_str}_R1"} = 1;
      
      if($is_paired_end){
        my $fastq_2 = "fastq_${pick_str}_2";
        $config->{$fastq_2} = {     
          "class" => "CQS::FilePickTask",
          "source_ref" => [$group_i],
          "sample_index" => 1, 
        };
        $config->{$task}{input_parameters}{"atac.fastqs_rep${pick_str}_R2_ref"} = [$fastq_2];
        $input_parameters_is_vector->{"atac.fastqs_rep${pick_str}_R2"} = 1;
      }

      $pick_index = $pick_index + 1;
    }

    $config->{$task}{input_parameters_is_vector} = $input_parameters_is_vector;
  }else{
    $config->{$task}{source_ref} = $files_ref;
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
    $config->{$task}{input_parameters}{"atac.fastqs_rep1_R1_ref"} = [$fastq_1];
    $config->{$task}{input_parameters}{"atac.fastqs_rep1_R2_ref"} = [$fastq_2];
  }

  push @$individual, $task;

  my $croo_task = $task . "_croo";
  $config->{$croo_task} = {
    class => "CQS::ProgramWrapperOneToOne",
    target_dir => "${target_dir}/${task_folder}_croo",
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
