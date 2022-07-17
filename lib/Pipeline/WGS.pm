#!/usr/bin/perl
package Pipeline::WGS;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use Pipeline::WdlPipeline;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performWGS 
  performWGSTask 
  add_bam_to_gvcf 
  add_gvcf_to_genotype
  add_gvcf_to_genotype_scatter 
  add_hard_filter_and_left_trim
  add_merge
  add_hard_filter_and_merge)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  #check files
  getValue($def, "ref_fasta");
  getValue($def, "dbsnp");
  getValue($def, "interval_list_file");
  getValue($def, "known_indels_sites_VCFs");

  #start from FASTQ
  initDefaultValue( $def, "perform_preprocessing", 1 );

  #let's start from bam file
  #initDefaultValue( $def, "perform_preprocessing", 0 );
  #initDefaultValue( $def, "bam_file_section", "files" );

  initDefaultValue( $def, "gatk_prefix", "gatk4_" );

  initDefaultValue( $def, "perform_replace_read_group", 0);
  initDefaultValue( $def, "replace_read_group_by_gatk4", 0);
  
  initDefaultValue( $def, "perform_mark_duplicates", 1);

  initDefaultValue( $def, "perform_gvcf_to_genotype", 1);

  initDefaultValue( $def, "perform_filter_and_merge", 1);
  initDefaultValue( $def, "use_hard_filter", "1" );
      
  initDefaultValue( $def, "max_thread", 8 );
  initDefaultValue( $def, "subdir",     0 );

  initDefaultValue( $def, "sra_to_fastq", 0 );

  initDefaultValue( $def, "aligner", "bwa" );
  if ( $def->{aligner} eq "bwa" ) {
    initDefaultValue( $def, "bwa_option", "" );
  }

  if(defined $def->{"ref_fasta_dict"} && (! defined $def->{"chromosome_names"})){
    my $dictFile = getValue($def, "ref_fasta_dict");
    my $primary_chromosome_only = getValue($def, "primary_chromosome_only", 1);
    my $chromosomes = readChromosomeFromDictFile($dictFile, $primary_chromosome_only);
    initDefaultValue( $def, "chromosome_names" , join(",", @$chromosomes));
  }

  initDefaultValue( $def, "filter_variants_by_allele_frequency",            1 );
  initDefaultValue( $def, "filter_variants_by_allele_frequency_percentage", 0.9 );
  initDefaultValue( $def, "filter_variants_by_allele_frequency_maf",        0.3 );
  initDefaultValue( $def, "filter_variants_fq_equal_1",                     0 );

  initDefaultValue( $def, "perform_cnv_gatk4_cohort", 1 );
  initDefaultValue( $def, "gatk4_cnv_by_scatter",     1 );
  initDefaultValue( $def, "gatk4_cnv_scatter_count",  100 );
  initDefaultValue( $def, "perform_cnv_gatk4_somatic", 0 );

  initDefaultValue( $def, "perform_muTect",       0 );
  initDefaultValue( $def, "perform_muTect2indel", 0 );
  initDefaultValue( $def, "perform_annovar",      0 );
  initDefaultValue( $def, "perform_cnv",          1 );
  initDefaultValue( $def, "perform_vep",          0 );
  initDefaultValue( $def, "perform_IBS",          0 );

  if ( $def->{perform_muTect} || $def->{perform_muTect2indel} ) {
    if ( defined $def->{mills} ) {
      initDefaultValue( $def, "indel_realignment", 1 );
      initDefaultValue( $def, "indel_vcf_files",   $def->{mills} );
    }
    else {
      initDefaultValue( $def, "indel_realignment", 0 );
    }
  }

  initDefaultValue( $def, "remove_duplicate", 1 );
  initDefaultValue( $def, "perform_multiqc", 0 );

  my $default_onco_options = {
    "picture_width" => "0",
    "picture_height" => "0",
    "sampleNamePattern" => ".",
    "BACKGROUND_color" => "lightgray",
    "MISSENSE_color" => "darkgreen",
    "MISSENSE_height" => 0.75,
    "TRUNC_color" => "red",
    "TRUNC_height" => 0.5, 
    "DUP_color" => "brown", 
    "DUP_height" => 0.25,
    "DEL_color" => "blue",
    "DEL_height" => 0.25,
    "MANUAL_order" => 0,
  };

  if (defined $def->{onco_options}) {
    $def->{onco_options} = merge_hash_left_precedent($def->{onco_options}, $default_onco_options);
  }else{
    $def->{onco_options} = $default_onco_options;
  }

  return $def;
}

sub add_bam_to_gvcf {
  my ($config, $def, $tasks, $target_dir, $gatk_prefix, $gatk_index_snv, $source ) = @_;

  my $bam_to_gvcf = {
    BaseRecalibratorScatter => {
      class             => "GATK4::BaseRecalibratorScatter",
      perform           => 1,
      target_dir        => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_BaseRecalibratorScatter",
      option            => "",
      source_ref        => $source,
      fasta_file        => $def->{ref_fasta},
      dbsnp_vcf         => $def->{dbsnp},
      known_indels_sites_VCFs => getValue($def, "known_indels_sites_VCFs"),
      sh_direct         => 0,
      pbs               => {
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    },
    GatherBQSRReports => {
      class                 => "GATK4::GatherBQSRReports",
      perform               => 1,
      target_dir            => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_GatherBQSRReports",
      option                => "",
      gather_name_ref       => $source,
      source_ref            => "BaseRecalibratorScatter",
      sh_direct             => 0,
      pbs                   => {
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "5gb"
      },
    },
    ApplyBQSRScatter => {
      class             => "GATK4::ApplyBQSRScatter",
      perform           => 1,
      target_dir        => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_ApplyBQSRScatter",
      option            => "",
      source_ref        => $source,
      fasta_file        => $def->{ref_fasta},
      bqsr_report_files_ref => "GatherBQSRReports",
      sh_direct         => 0,
      pbs               => {
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    },
    GatherSortedBamFiles => {
      class                 => "GATK4::GatherSortedBamFilesSambamba",
      perform               => 1,
      target_dir            => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_GatherSortedBamFiles",
      option                => "",
      gather_name_ref       => $source,
      source_ref            => ["ApplyBQSRScatter"],
      extension             => ".recalibrated.bam",
      sh_direct             => 0,
      pbs                   => {
        "nodes"    => "1:ppn=8",
        "walltime" => "24",
        "mem"      => "80gb"
      },
    },
    HaplotypeCallerScatter => {
      class             => "GATK4::HaplotypeCallerScatter",
      perform           => 1,
      target_dir        => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_HaplotypeCallerScatter",
      #https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/haplotypecaller-gvcf-gatk4.wdl
      option            => "\\
    -contamination 0 \\
    -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \\
    -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90",
      source_ref        => ["GatherSortedBamFiles", ".bam\$"],
      java_option       => "",
      fasta_file        => $def->{ref_fasta},
      extension         => ".g.vcf.gz",
      blacklist_file    => $def->{blacklist_file},
      gvcf              => 1,
      sh_direct         => 0,
      pbs               => {
        "nodes"    => "1:ppn=8",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    },
    GatherVcfs => {
      class                 => "GATK4::GatherVcfs",
      perform               => 1,
      target_dir            => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_GatherVcfs",
      option                => "",
      gather_name_ref       => $source,
      source_ref            => ["HaplotypeCallerScatter"],
      extension             => ".g.vcf.gz",
      sh_direct             => 0,
      pbs                   => {
        "nodes"    => "1:ppn=1",
        "walltime" => "4",
        "mem"      => "10gb"
      },
    },
  };

  my @newtasks = ("BaseRecalibratorScatter", 
      "GatherBQSRReports", 
      "ApplyBQSRScatter", 
      "GatherSortedBamFiles",
      "HaplotypeCallerScatter", 
      "GatherVcfs");

  push(@$tasks, @newtasks);

  foreach my $task (@newtasks){
    $config->{$task} = $bam_to_gvcf->{$task};
  }

  return "GatherVcfs";
}

sub add_gvcf_to_genotype_scatter {
  my ($config, $def, $tasks, $target_dir, $gatk_prefix, $gatk_index_snv, $source ) = @_;

  my $gvcf_to_genotype = {
    GenomicsDBImportScatter => {
      class             => "GATK4::GenomicsDBImportScatter",
      perform           => 1,
      target_dir        => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_GenomicsDBImportScatter",
      option            => "",
      source_ref        => $source,
      java_option       => "",
      sh_direct         => 0,
      pbs               => {
        "nodes"    => "1:ppn=4",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    },
    GenotypeGVCFsScatter => {
      class             => "GATK4::GenotypeGVCFsScatter",
      perform           => 1,
      target_dir        => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_GenotypeGVCFsScatter",
      option            => "",
      source_ref        => ["GenomicsDBImportScatter"],
      fasta_file        => $def->{ref_fasta},
      dbsnp_vcf         => $def->{dbsnp},
      java_option       => "",
      sh_direct         => 0,
      pbs               => {
        "nodes"    => "1:ppn=4",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    }
  };

  my @newtasks = ("GenomicsDBImportScatter", "GenotypeGVCFsScatter");

  push(@$tasks, @newtasks);

  foreach my $task (@newtasks){
    $config->{$task} = $gvcf_to_genotype->{$task};
  }

  return "GenotypeGVCFsScatter";
}

sub add_gvcf_to_genotype {
  my ($config, $def, $tasks, $target_dir, $gatk_prefix, $gatk_index_snv, $source, $interval_key ) = @_;

  my $gvcf_to_genotype = {
    GenomicsDBImport => {
      class             => "GATK4::GenomicsDBImport",
      perform           => 1,
      target_dir        => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_GenomicsDBImport",
      option            => "",
      source_ref        => $source,
      target_intervals_file => getValue($def, $interval_key),
      java_option       => "",
      sh_direct         => 0,
      pbs               => {
        "nodes"    => "1:ppn=4",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    },
    GenotypeGVCFs => {
      class             => "GATK4::GenotypeGVCFs",
      perform           => 1,
      target_dir        => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_GenotypeGVCFs",
      option            => "",
      source_ref        => ["GenomicsDBImport"],
      fasta_file        => $def->{ref_fasta},
      dbsnp_vcf         => $def->{dbsnp},
      target_intervals_file => getValue($def, $interval_key),
      java_option       => "",
      sh_direct         => 0,
      pbs               => {
        "nodes"    => "1:ppn=4",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    }
  };

  my @newtasks = ("GenomicsDBImport", "GenotypeGVCFs");

  push(@$tasks, @newtasks);

  foreach my $task (@newtasks){
    $config->{$task} = $gvcf_to_genotype->{$task};
  }

  return "GenotypeGVCFs";
}

sub add_hard_filter_and_left_trim {
  my ($config, $def, $tasks, $target_dir, $gatk_prefix, $gatk_index_snv, $source ) = @_;

  my $filter = {
    VariantFilterHard  => {
      class                 => "GATK4::VariantFilterHard",
      perform               => 1,
      target_dir            => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_VariantFilterHard",
      source_ref            => $source,
      option                => "",
      docker_prefix         => "gatk4_",
      is_sample_size_small  => 1,
      sh_direct             => 1,
      pbs                   => {
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    },
    LeftTrim  => {
      class                 => "GATK4::LeftTrim",
      perform               => 1,
      target_dir            => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_LeftTrim",
      source_ref            => "VariantFilterHard",
      option                => "",
      docker_prefix         => "cqs_",
      fasta_file            => $def->{ref_fasta},
      extension             => ".indels.snp.hardfilter.pass.norm.nospan.vcf.gz",
      sh_direct             => 0,
      pbs                   => {
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      },
    },
  };

  my @newtasks = ("VariantFilterHard", "LeftTrim");

  push(@$tasks, @newtasks);

  foreach my $task (@newtasks){
    $config->{$task} = $filter->{$task};
  }

  return "LeftTrim";
}

sub add_merge {
  my ($config, $def, $tasks, $target_dir, $gatk_prefix, $gatk_index_snv, $source ) = @_;

  my $task_name = "MergeVcfs";
  $config->{$task_name} = {
    class                 => "GATK4::MergeVcfs",
    perform               => 1,
    target_dir            => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_MergeVcfs",
    option                => "",
    source_ref            => ["LeftTrim"],
    extension             => ".vcf.gz",
    docker_prefix         => "gatk4_",
    sh_direct             => 0,
    pbs                   => {
      "nodes"    => "1:ppn=1",
      "walltime" => "4",
      "mem"      => "10gb"
    },
  };

  push(@$tasks, $task_name);

  return $task_name;
}

sub add_hard_filter_and_merge {
  my ($config, $def, $tasks, $target_dir, $gatk_prefix, $gatk_index_snv, $source ) = @_;

  my $filter = add_hard_filter_and_left_trim($config, $def, $tasks, $target_dir, $gatk_prefix, $gatk_index_snv, $source);
  my $merge = add_merge($config, $def, $tasks, $target_dir, $gatk_prefix, $gatk_index_snv, $filter);

  return $merge;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  $config->{general}{interval_list_file} = getValue($def, "interval_list_file");

  my $gatk_prefix = getValue($def, "gatk_prefix");
  my $gatk_index_snv = "SNV_index";
  my $use_hard_filter = getValue($def, "use_hard_filter");

  my $bam_section = $def->{bam_file_section};
  if (not defined $bam_section) {
    my $aligner_scatter_count = getValue($def, "aligner_scatter_count", 0);

    my $bam_ref;
    my $bam_task;
    if ($aligner_scatter_count > 0 || getValue($def, "has_multiple_fastq_per_sample", 0)){
      my ($bam_ref, $bam_task) = add_BWA_and_summary_scatter($config, $def, $individual, $target_dir, "files");
      $bam_section = $bam_ref;
    } else {
      my $bam_task = "bwa_wgs";
      add_BWA_WGS($config, $def, $individual, $target_dir, $bam_task, "files");
      $bam_section = [ $bam_task, ".bam\$" ];
    }
  }

  my $perform_replace_read_group = getValue($def, "perform_replace_read_group", 0);
  if ($perform_replace_read_group) {

    if ($def->{replace_read_group_by_gatk4}){
      $config->{replace_read_group} = {
        class                 => "CQS::ProgramWrapperOneToOne",
        perform               => 1,
        target_dir            => "${target_dir}/". $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_replace_read_group",
        option                => "AddOrReplaceReadGroups --CREATE_INDEX true --TMP_DIR tmp -I __FILE__ -O __NAME__.tmp.bam \\
  -ID 1 \\
  -LB __NAME__ \\
  -PL ILLUMINA \\
  -PU __NAME__ \\
  -SM __NAME__ && touch __OUTPUT__.done

if [[ -e __OUTPUT__.done ]]; then
  mv __NAME__.tmp.bam __OUTPUT__
  mv __NAME__.tmp.bai __OUTPUT__.bai
  rm __OUTPUT__.done
fi
",
        interpretor           => "",
        program               => "gatk",
        docker_prefix         => "gatk4_",
        check_program         => 0,
        source_arg            => "-i",
        source_ref            => $bam_section,
        source_join_delimiter => "",
        output_to_same_folder => 1,
        output_arg            => "-o",
        output_file_prefix    => ".rg.bam",
        output_file_ext       => ".rg.bam",
        sh_direct             => 0,
        pbs                   => {
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "40gb"
        },
      };
    }else{
      $config->{replace_read_group} = {
        class                 => "CQS::ProgramWrapperOneToOne",
        perform               => 1,
        target_dir            => "${target_dir}/". $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_replace_read_group",
        option                => getValue($def, "replace_read_group_option"),
        interpretor           => "python3",
        program               => "../Format/replace_rg.py",
        docker_prefix         => "cqs_",
        check_program         => 1,
        source_arg            => "-i",
        source_ref            => $bam_section,
        source_join_delimiter => "",
        output_to_same_folder => 1,
        output_arg            => "-o",
        output_file_prefix    => ".rg.bam",
        output_file_ext       => ".rg.bam",
        sh_direct             => 0,
        pbs                   => {
          "nodes"    => "1:ppn=8",
          "walltime" => "24",
          "mem"      => "40gb"
        },
      };
    }
    $bam_section = "replace_read_group";
    push(@$summary, $bam_section);
  }

  my $perform_mark_duplicates = getValue($def, "perform_mark_duplicates", 1);
  if ($perform_mark_duplicates) {
    my $markduplicates = "markduplicates";

    addMarkduplicates($config, $def, $summary, $target_dir, $markduplicates, $bam_section);

    $bam_section = "markduplicates";
  }

  my $gvcf_section = add_bam_to_gvcf($config, $def, $summary, $target_dir, $gatk_prefix, $gatk_index_snv, $bam_section);

  if ($def->{perform_gvcf_to_genotype}) {
    my $genotypeGVCFs_section = add_gvcf_to_genotype_scatter($config, $def, $summary, $target_dir, $gatk_prefix, $gatk_index_snv, $gvcf_section);

    if ($def->{perform_filter_and_merge}) {
      my $merged_vcf_section;
      if (not $use_hard_filter){
    # We will have to do hard filter directly for mouse dataset
    #   #https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/JointGenotyping.wdl
    #   HardFilterAndMakeSitesOnlyVcf => {
    #     class                 => "CQS::ProgramWrapperOneToOne",
    #     perform               => 1,
    #     target_dir            => "${target_dir}/nih_bam_10_HardFilterAndMakeSitesOnlyVcf",
    #     option                => "--java-options \"-Xms10g -Xmx10g\" \\
    #   VariantFiltration \\
    #   --filter-expression \"ExcessHet > 54.69\" \\
    #   --filter-name ExcessHet \\
    #   -O __NAME__.variant_filtered.vcf.gz \\
    #   -V __FILE__

    # gatk --java-options \"-Xms10g -Xmx10g\" \\
    #   MakeSitesOnlyVcf \\
    #   -I __NAME__.variant_filtered.vcf.gz \\
    #   -O __NAME__.sites_only.variant_filtered.vcf.gz

    # ",
    #     interpretor           => "",
    #     program               => "gatk",
    #     docker_prefix         => "gatk4_",
    #     check_program         => 0,
    #     source_arg            => "",
    #     source_ref            => ["GenotypeGVCFs"],
    #     source_join_delimiter => " ",
    #     output_to_same_folder => 1,
    #     output_arg            => "-O",
    #     output_file_prefix    => "",
    #     output_file_ext       => ".variant_filtered.vcf.gz",
    #     output_other_ext      => ".sites_only.variant_filtered.vcf.gz",
    #     sh_direct             => 0,
    #     pbs                   => {
    #       "nodes"    => "1:ppn=1",
    #       "walltime" => "24",
    #       "mem"      => "20gb"
    #     },
    #   },
    #   SitesOnlyGatherVcf => {
    #     class                 => "CQS::ProgramWrapper",
    #     perform               => 1,
    #     target_dir            => "${target_dir}/nih_bam_11_SitesOnlyGatherVcf",
    #     option                => "--java-options -Xms6g GatherVcfsCloud  \\
    #   --input __FILE__ \\
    #   --output __NAME__.sites_only.vcf.gz && tabix __NAME__.sites_only.vcf.gz",
    #     interpretor           => "",
    #     program               => "gatk",
    #     docker_prefix         => "gatk4_",
    #     check_program         => 0,
    #     source_arg            => "--input",
    #     source_type           => "array",
    #     source_ref            => ["HardFilterAndMakeSitesOnlyVcf"],
    #     source_join_delimiter => " \\\n  --input ",
    #     output_arg            => "--output",
    #     output_file_prefix    => "",
    #     output_file_ext       => ".sites_only.vcf.gz",
    #     sh_direct             => 0,
    #     pbs                   => {
    #       "nodes"    => "1:ppn=1",
    #       "walltime" => "4",
    #       "mem"      => "10gb"
    #     },
    #   },
      } else {
        $merged_vcf_section = add_hard_filter_and_merge($config, $def, $summary, $target_dir, $gatk_prefix, $gatk_index_snv, $genotypeGVCFs_section);
      }

      addAnnovar($config, $def, $summary, $target_dir, $merged_vcf_section, , "", $gatk_prefix, $def, $gatk_index_snv );
    }
  }

  my $tasks = [@$individual, @$summary];

  $config->{sequencetask} = {
    class => "CQS::SequenceTaskSlurmSlim",
    perform => 1,
    target_dir            => "${target_dir}/sequencetask",
    option                => "",
    source                => {
      
      processing => $tasks,
    },
    sh_direct             => 1,
    pbs                   => {
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  };

  return($config);
};

sub performWGS {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

sub performWGSTask {
  my ( $def, $task ) = @_;
  my $config = getConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
