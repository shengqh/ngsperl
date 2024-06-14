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

our %EXPORT_TAGS = ( 'all' => [qw(
  initializeWGSDefaultOptions
  performWGS 
  performWGSTask 
  add_bam_recalibration
  add_recalibrated_bam_to_gvcf
  add_bam_to_gvcf 
  add_gvcf_to_genotype
  add_gvcf_to_genotype_scatter 
  add_hard_filter_and_left_trim
  add_merge
  add_hard_filter_and_merge)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeWGSDefaultOptions {
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

  initDefaultValue( $def, "callvariants_vqsr_mode", 1 );
  
  initDefaultValue( $def, "max_thread", 8 );
  initDefaultValue( $def, "subdir",     0 );

  initDefaultValue( $def, "sra_to_fastq", 0 );

  initDefaultValue( $def, "aligner", "bwa" );
  if ( $def->{aligner} eq "bwa" ) {
    initDefaultValue( $def, "bwa_option", "-K 100000000 -v 3" );
  }

  initDefaultValue( $def, "mark_duplicates_use_tmp_folder", 0 );

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

  initDefaultValue( $def, "recalibration_annotation_values", ["QD", 
  "MQRankSum",
  "ReadPosRankSum",
  "FS",
  "MQ",
  "SOR",
  "DP"]);

  initDefaultValue( $def, "snp_recalibration_tranche_values", ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0" ]);
  initDefaultValue( $def, "snp_recalibration_annotation_values", ["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"]);
  initDefaultValue( $def, "indel_recalibration_tranche_values", ["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"]);
  initDefaultValue( $def, "indel_recalibration_annotation_values",  ["FS", "ReadPosRankSum", "MQRankSum", "QD", "SOR", "DP"]);
  initDefaultValue( $def, "snp_filter_level", 99.7 );
  initDefaultValue( $def, "indel_filter_level", 99.7 );
  initDefaultValue( $def, "SNP_VQSR_downsampleFactor", 10 );

  initDefaultValue( $def, "perform_muTect",       0 );
  initDefaultValue( $def, "perform_muTect2indel", 0 );

  initDefaultValue( $def, "perform_annovar",      1 );
  initDefaultValue( $def, "maximum_freq_values",  "0.05,0.01" );

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

  my $perform_split_fastq = getValue($def, "perform_split_fastq", 0);
  if($perform_split_fastq) {#[ 0, "by_dynamic", "by_file", "by_scatter"],
    if ($perform_split_fastq eq "by_dynamic"){ 
      initDefaultValue($def, "split_fastq_min_file_size_gb", 10); #only the file with file size larger than this number would be splitted
      initDefaultValue($def, "split_fastq_trunk_file_size_gb", 5); #the splitted file will be smaller than this file size
    }elsif($perform_split_fastq eq "by_scatter"){ 
      initDefaultValue($def, "aligner_scatter_count", 10); #split data into equal number of small files
    }elsif($perform_split_fastq eq "by_file"){ 
      #nothing need to be set
    }else{
      die 'wrong perform_split_fastq value, it should be one of [ 0, "by_dynamic", "by_file", "by_scatter"]';
    }
  }

  return $def;
}

sub add_bam_recalibration {
  my ($config, $def, $tasks, $target_dir, $gatk_prefix, $gatk_index_snv, $source ) = @_;

  my $to_cram = getValue($def, "GatherSortedBamFiles_cram", 0);

  my $bam_suffix = $to_cram ? ".cram": ".bam";
  my $bam_index_suffix = $to_cram ? ".cram.crai": ".bai";

  my $bam_recalibration = {
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
      class                 => "GATK4::GatherSortedBamFiles",
      perform               => 1,
      target_dir            => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_GatherSortedBamFiles",
      option                => "",
      gather_name_ref       => $source,
      source_ref            => ["ApplyBQSRScatter"],
      ref_fasta => getValue($def, "ref_fasta"),
      sh_direct             => 0,
      pbs                   => {
        "nodes"    => "1:ppn=8",
        "walltime" => "24",
        "mem"      => "80gb"
      },
    },
  };
  my @newtasks = ("BaseRecalibratorScatter", 
      "GatherBQSRReports", 
      "ApplyBQSRScatter", 
      "GatherSortedBamFiles");

  push(@$tasks, @newtasks);
  my $bam_task = "GatherSortedBamFiles";

  if($to_cram){
    my $BamToCram = "BamToCram";
    my $ref_fasta = getValue($def, "ref_fasta");
    $config->{$BamToCram} = {
      class              => "CQS::ProgramWrapperOneToOne",
      perform            => 1,
      target_dir => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_${BamToCram}",
      program => "",
      check_program => 0,
      source_ref => "GatherSortedBamFiles",
      option => "
echo sort_bam_to_cram=`date`
samtools sort -m 5G \\
  --output-fmt CRAM \\
  --reference $ref_fasta \\
  --threads 8 \\
  --write-index \\
  -o tmp.__OUTPUT__ __FILE__

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.succeed 
  rm -f __NAME__.failed 

  mv tmp.__OUTPUT__ __OUTPUT__
  mv tmp.__OUTPUT__.crai __OUTPUT__.crai
else
  touch __NAME__.failed
  rm -f __NAME__.succeed tmp.__OUTPUT__ tmp.__OUTPUT__.crai
fi
",
      output_file_prefix => ".cram",
      output_file_ext => ".cram",
      output_by_file => 0,
      use_tmp_folder => 0,
      sh_direct          => 0,
      pbs                => {
        "nodes"     => "1:ppn=8",
        "walltime"  => "24",
        "mem"       => "40gb"
      },
    };
    push(@$tasks, $BamToCram);

    $bam_task = $BamToCram;
  }

  foreach my $task (@newtasks){
    $config->{$task} = $bam_recalibration->{$task};
  }

  return $bam_task;
}

sub add_recalibrated_bam_to_gvcf {
  my ($config, $def, $tasks, $target_dir, $gatk_prefix, $gatk_index_snv, $source ) = @_;

  my $bam_to_gvcf = {
    HaplotypeCallerScatter => {
      class             => "GATK4::HaplotypeCallerScatter",
      perform           => 1,
      target_dir        => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_HaplotypeCallerScatter",
      #https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/haplotypecaller-gvcf-gatk4.wdl
      option            => "\\
  -contamination 0 \\
  -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \\
  -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \\
  --soft-clip-low-quality-ends true \\
  --dont-use-soft-clipped-bases true",
      source_ref        => $source,
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

  my @newtasks = ("HaplotypeCallerScatter", "GatherVcfs");

  push(@$tasks, @newtasks);

  foreach my $task (@newtasks){
    $config->{$task} = $bam_to_gvcf->{$task};
  }

  return "GatherVcfs";
}

sub add_bam_to_gvcf {
  my ($config, $def, $tasks, $target_dir, $gatk_prefix, $gatk_index_snv, $source ) = @_;
  my $bam_recalibration_section = add_bam_recalibration($config, $def, $tasks, $target_dir, $gatk_prefix, $gatk_index_snv, $source);
  my $gvcf_section = add_recalibrated_bam_to_gvcf($config, $def, $tasks, $target_dir, $gatk_prefix, $gatk_index_snv, $bam_recalibration_section);
  return($gvcf_section);
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
        "nodes"    => "1:ppn=6",
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
        "walltime" => "4",
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

  $def = initializeWGSDefaultOptions($def);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);
  #merge summary and individual 
  push @$individual, @$summary;
  $summary = $individual;

  if(defined $def->{annotation_genes}){
    addGeneLocus($config, $def, $summary, $target_dir);
  }

  $config->{general}{interval_list_file} = getValue($def, "interval_list_file");

  my $gatk_prefix = getValue($def, "gatk_prefix");
  my $gatk_index_snv = "SNV_index";
  my $callvariants_vqsr_mode = getValue($def, "callvariants_vqsr_mode");

  my $bam_section = undef;
  if (getValue($def, "is_input_bam", 0)){
    $bam_section = $def->{files};
  }elsif(defined $def->{bam_file_section}){
    $bam_section = $def->{bam_file_section};
  }else{
    my $perform_split_fastq = getValue($def, "perform_split_fastq", 0);

    my $bam_ref;
    my $bam_task;
    if ($perform_split_fastq){
      my ($bam_ref, $bam_task) = add_BWA_and_summary_scatter($config, $def, $individual, $target_dir, $source_ref);
      $bam_section = $bam_ref;
    } else {
      my $bam_task = "bwa_wgs";
      add_BWA_WGS($config, $def, $individual, $target_dir, $bam_task, $source_ref);
      $bam_section = [ $bam_task, ".bam\$" ];
    }
  }

  my $ref_fasta = getValue($def, "ref_fasta");

  my $merged_vcf_section = undef;
  if (getValue($def, "is_input_vcf", 0)){
    $merged_vcf_section = "files";

    if(getValue($def, "filter_input_vcf", 1)){
      $config->{filter_input_vcf} = {
        class                 => "CQS::ProgramWrapperOneToOne",
        perform               => 1,
        target_dir            => "${target_dir}/". $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_filter_input_vcf",
        option                => "
rm -f __OUTPUT__.failed

gatk --java-options '-Xmx60g' SelectVariants \\
  -R $ref_fasta \\
  -V __FILE__ \\
  -O __OUTPUT__ --exclude-filtered

status=\$?
if [[ \$status -ne 0 ]]; then
  touch __OUTPUT__.failed
  rm -f __OUTPUT__ __OUTPUT__.tbi
else
  touch __OUTPUT__.succeed
fi  
",
        interpretor           => "",
        program               => "",
        docker_prefix         => "gatk4_",
        check_program         => 0,
        source_arg            => "-i",
        source_ref            => $merged_vcf_section,
        source_join_delimiter => "",
        output_to_same_folder => 1,
        output_arg            => "-o",
        output_file_prefix    => ".pass.vcf.gz",
        output_file_ext       => ".pass.vcf.gz,.pass.vcf.gz.tbi",
        sh_direct             => 0,
        pbs                   => {
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "60gb"
        },
      };
      $config->{LeftTrim} = {
        class                 => "GATK4::LeftTrim",
        perform               => 1,
        target_dir            => "${target_dir}/" . $gatk_prefix . getNextIndex($def, $gatk_index_snv) . "_LeftTrim",
        source_ref            => "filter_input_vcf",
        option                => "",
        docker_prefix         => "cqs_",
        fasta_file            => $ref_fasta,
        extension             => ".pass.norm.nospan.vcf.gz",
        sh_direct             => 0,
        pbs                   => {
          "nodes"    => "1:ppn=1",
          "walltime" => "10",
          "mem"      => "40gb"
        },
      };
      $merged_vcf_section = "LeftTrim";
    }
  }else{
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
      if(!defined $config->{bwa_02_markduplicates}){
        my $perform_split_fastq = getValue($def, "perform_split_fastq", 0);
        #if perform_split_fastq, markduplicates might be added with split_fastq tasks.
        my $markduplicates = "markduplicates";
        addMarkduplicates($config, $def, $summary, $target_dir, $markduplicates, $bam_section);
        $bam_section = "markduplicates";
      }
    }

    my $bam_recalibration_section = add_bam_recalibration($config, $def, $summary, $target_dir, $gatk_prefix, $gatk_index_snv, $bam_section);
    
    my $gvcf_section = add_recalibrated_bam_to_gvcf($config, $def, $summary, $target_dir, $gatk_prefix, $gatk_index_snv, $bam_recalibration_section);

    if($def->{perform_extract_bam}){
      my $extract_bam_locus = getValue($def, "extract_bam_locus");
      my $extract_bam_task = "extract_bam_locus" . getValue($def, "extract_bam_locus_suffix", "");
      add_extract_bam_locus($config, $def, $individual, $target_dir, $extract_bam_task, $extract_bam_locus, $bam_recalibration_section );
    }

    #my $gvcf_section = add_bam_to_gvcf($config, $def, $summary, $target_dir, $gatk_prefix, $gatk_index_snv, $bam_section);

    if ( $def->{perform_cnv_gatk4_cohort} ) {
      my $cnv_prefix = getValue($def, "cnv_prefix", "");
      my $cnvMap = addGATK4CNVGermlineCohortAnalysis( $config, $def, $target_dir, $bam_recalibration_section, $cnv_prefix, $individual, $summary, $summary, $summary, $summary, $summary );
    }

    if ($def->{perform_gvcf_to_genotype}) {
      my $genotypeGVCFs_section = add_gvcf_to_genotype_scatter($config, $def, $summary, $target_dir, $gatk_prefix, $gatk_index_snv, $gvcf_section);

      if ($def->{perform_filter_and_merge}) {
        if ($callvariants_vqsr_mode){
          #stop("not implemented");
          
          my $HardFilterAndMakeSitesOnlyVcf = "HardFilterAndMakeSitesOnlyVcf";
          $config->{$HardFilterAndMakeSitesOnlyVcf} = {
            class                 => "CQS::ProgramWrapperOneToOne",
            perform               => 1,
            target_dir            => "${target_dir}/" . $gatk_prefix .  getNextIndex($def, $gatk_index_snv) . "_HardFilterAndMakeSitesOnlyVcf",
            option                => "
echo VariantFiltration ...
gatk --java-options \"-Xms10g -Xmx18g\" \\
  VariantFiltration \\
  --filter-expression \"ExcessHet > 54.69\" \\
  --filter-name ExcessHet \\
  -O __NAME__.variant_filtered.vcf.gz \\
  -V __FILE__

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.variant_filtered.vcf.gz.succeed
  rm -f __NAME__.variant_filtered.vcf.gz.failed

  echo MakeSitesOnlyVcf ...
  gatk --java-options \"-Xms10g -Xmx18g\" \\
    MakeSitesOnlyVcf \\
    -I __NAME__.variant_filtered.vcf.gz \\
    -O __NAME__.tmp.variant_filtered.sites_only.vcf.gz

  status=\$?
  if [[ \$status -eq 0 ]]; then
    touch __NAME__.variant_filtered.sites_only.vcf.gz.succeed
    rm -f __NAME__.variant_filtered.sites_only.vcf.gz.failed
    mv __NAME__.tmp.variant_filtered.sites_only.vcf.gz __NAME__.variant_filtered.sites_only.vcf.gz
    mv __NAME__.tmp.variant_filtered.sites_only.vcf.gz.tbi __NAME__.variant_filtered.sites_only.vcf.gz.tbi
  else
    touch __NAME__.variant_filtered.sites_only.vcf.gz.failed
    rm -f __NAME__.variant_filtered.sites_only.vcf.gz.succeed __NAME__.tmp.variant_filtered.sites_only.vcf.gz __NAME__.tmp.variant_filtered.sites_only.vcf.gz.tbi
  fi
else
  touch __NAME__.variant_filtered.vcf.gz.failed
  rm -f __NAME__.variant_filtered.vcf.gz.succeed __NAME__.variant_filtered.vcf.gz __NAME__.variant_filtered.vcf.gz.tbi
fi
",
            interpretor           => "",
            program               => "",
            docker_prefix         => "gatk4_",
            check_program         => 0,
            source_arg            => "",
            source_ref            => $genotypeGVCFs_section,
            other_localization_ext_array => [".tbi"],
            source_join_delimiter => " ",
            output_to_same_folder => 1,
            output_arg            => "-O",
            output_file_prefix    => "",
            output_file_ext       => ".variant_filtered.vcf.gz",
            output_other_ext      => ".variant_filtered.sites_only.vcf.gz",
            sh_direct             => 0,
            pbs                   => {
              "nodes"    => "1:ppn=1",
              "walltime" => "24",
              "mem"      => "20gb"
            },
          };

          my $SitesOnlyGatherVcf = "SitesOnlyGatherVcf";
          $config->{$SitesOnlyGatherVcf} = {
            class                 => "CQS::ProgramWrapper",
            perform               => 1,
            target_dir            => "${target_dir}/" . $gatk_prefix .  getNextIndex($def, $gatk_index_snv) . "_SitesOnlyGatherVcf",
            option                => "
# --ignore-safety-checks makes a big performance difference so we include it in our invocation.
# This argument disables expensive checks that the file headers contain the same set of
# genotyped samples and that files are in order by position of first record.
gatk --java-options \"-Xms6000m -Xmx6500m\" \\
  GatherVcfsCloud \\
  --ignore-safety-checks \\
  --gather-type BLOCK \\
  --input __FILE__ \\
  --output __NAME__.tmp.sites_only.vcf.gz 
  
status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.sites_only.vcf.gz.succeed
  rm -f __NAME__.sites_only.vcf.gz.failed
  mv __NAME__.tmp.sites_only.vcf.gz __NAME__.sites_only.vcf.gz 
  mv __NAME__.tmp.sites_only.vcf.gz.tbi __NAME__.sites_only.vcf.gz.tbi

  tabix __NAME__.sites_only.vcf.gz
else
  touch __NAME__.sites_only.vcf.gz.failed
  rm -f __NAME__.sites_only.vcf.gz.succeed __NAME__.tmp.sites_only.vcf.gz __NAME__.tmp.sites_only.vcf.gz.tbi
fi
",
            interpretor           => "",
            program               => "",
            docker_prefix         => "gatk4_",
            check_program         => 0,
            source_arg            => "--input",
            source_type           => "array",
            source_ref            => [ $HardFilterAndMakeSitesOnlyVcf, ".sites_only.vcf.gz" ],
            source_join_delimiter => " \\\n  --input ",
            output_arg            => "--output",
            output_file_prefix    => "",
            output_file_ext       => ".sites_only.vcf.gz",
            sh_direct             => 0,
            pbs                   => {
              "nodes"    => "1:ppn=1",
              "walltime" => "24",
              "mem"      => "7gb"
            },
          };

          my $mills_resource_vcf = getValue($def, "mills");
          my $dbsnp_resource_vcf = getValue($def, "dbsnp");
          my $axiomPoly_resource_vcf = getValue($def, "axiomPoly");

          my $indel_tranche = getValue($def, "indel_recalibration_tranche_values");
          my $indel_tranche_str = join(" -tranche ", @$indel_tranche);

          my $indel_anno = getValue($def, "indel_recalibration_annotation_values");
          my $indel_anno_str = join(" -an ", @$indel_anno);

          my $indel_max_gaussians = getValue($def, "indel_max_gaussians", 4);

          my $IndelsVariantRecalibrator = "IndelsVariantRecalibrator";
          $config->{$IndelsVariantRecalibrator} = {
            class                 => "CQS::ProgramWrapperOneToOne",
            perform               => 1,
            target_dir            => "${target_dir}/" . $gatk_prefix .  getNextIndex($def, $gatk_index_snv) . "_IndelsVariantRecalibrator",
            option                => "
gatk --java-options \"-Xms24g -Xms25g\" \\
  VariantRecalibrator \\
  -AS \\
  -V __FILE__ \\
  -O __NAME__.tmp.indels.recal.vcf.gz \\
  --tranches-file __NAME__.indels.tranches \\
  --trust-all-polymorphic \\
  -tranche $indel_tranche_str \\
  -an $indel_anno_str \\
  -mode INDEL \\
  --max-gaussians $indel_max_gaussians \\
  --resource:mills,known=false,training=true,truth=true,prior=12 ${mills_resource_vcf} \\
  --resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiomPoly_resource_vcf} \\
  --resource:dbsnp,known=true,training=false,truth=false,prior=2 ${dbsnp_resource_vcf}
 
status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.indels.recal.vcf.gz.succeed
  rm -f __NAME__.indels.recal.vcf.gz.failed
  mv __NAME__.tmp.indels.recal.vcf.gz __NAME__.indels.recal.vcf.gz
  mv __NAME__.tmp.indels.recal.vcf.gz.tbi __NAME__.indels.recal.vcf.gz.tbi
else
  touch __NAME__.indels.recal.vcf.gz.failed
  rm -f __NAME__.indels.recal.vcf.gz.succeed __NAME__.tmp.indels.recal.vcf.gz __NAME__.tmp.indels.recal.vcf.gz.tbi
fi
",
            interpretor           => "",
            program               => "",
            docker_prefix         => "gatk4_",
            check_program         => 0,
            source_arg            => "-V",
            source_ref            => $SitesOnlyGatherVcf,
            output_arg            => "-O",
            output_file_prefix    => "",
            output_file_ext       => ".indels.recal.vcf.gz",
            output_other_ext      => ".indels.tranches",
            output_to_same_folder => 1,
            no_output             => 1,
            other_localization_ext_array => [".tbi"],
            sh_direct             => 0,
            pbs                   => {
              "nodes"    => "1:ppn=2",
              "walltime" => "12",
              "mem"      => "26gb"
            },
          };

          my $hapmap_resource_vcf = getValue($def, "hapmap");
          my $omni_resource_vcf = getValue($def, "omni");
          my $one_thousand_genomes_resource_vcf = getValue($def, "g1000");

          my $snp_tranche = getValue($def, "snp_recalibration_tranche_values");
          my $snp_tranche_str = join(" -tranche ", @$snp_tranche);

          my $snp_anno = getValue($def, "snp_recalibration_annotation_values");
          my $snp_anno_str = join(" -an ", @$snp_anno);

          my $snp_max_gaussians = getValue($def, "snp_max_gaussians", 6);

          my $SNPsVariantRecalibrator = "SNPsVariantRecalibratorClassic";
          $config->{$SNPsVariantRecalibrator} = {
            class                 => "CQS::ProgramWrapperOneToOne",
            perform               => 1,
            target_dir            => "${target_dir}/" . $gatk_prefix .  getNextIndex($def, $gatk_index_snv) . "_SNPsVariantRecalibratorClassic",
            option                => "
gatk --java-options \"-Xmx24g -Xms24g\" \\
  VariantRecalibrator \\
  -AS \\
  -V __FILE__ \\
  -O __NAME__.tmp.snps.recal.vcf.gz \\
  --tranches-file __NAME__.snps.tranches \\
  --trust-all-polymorphic \\
  -tranche $snp_tranche_str \\
  -an $snp_anno_str \\
  -mode SNP \\
  --max-gaussians $snp_max_gaussians \\
  --resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap_resource_vcf} \\
  --resource:omni,known=false,training=true,truth=true,prior=12 ${omni_resource_vcf} \\
  --resource:1000G,known=false,training=true,truth=false,prior=10 ${one_thousand_genomes_resource_vcf} \\
  --resource:dbsnp,known=true,training=false,truth=false,prior=7 ${dbsnp_resource_vcf}

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.snps.recal.vcf.gz.succeed
  rm -f __NAME__.snps.recal.vcf.gz.failed
  mv __NAME__.tmp.snps.recal.vcf.gz __NAME__.snps.recal.vcf.gz
  mv __NAME__.tmp.snps.recal.vcf.gz.tbi __NAME__.snps.recal.vcf.gz.tbi
else
  touch __NAME__.snps.recal.vcf.gz.failed
  rm -f __NAME__.snps.recal.vcf.gz.succeed __NAME__.tmp.snps.recal.vcf.gz __NAME__.tmp.snps.recal.vcf.gz.tbi
fi
",
            interpretor           => "",
            program               => "",
            docker_prefix         => "gatk4_",
            check_program         => 0,
            source_arg            => "-V",
            source_ref            => $SitesOnlyGatherVcf,
            output_arg            => "-O",
            output_file_prefix    => "",
            output_file_ext => ".snps.recal.vcf.gz",
            output_other_ext      => ".snps.tranches",
            output_to_same_folder => 1,
            no_output => 1,
            other_localization_ext_array => [".tbi"],
            sh_direct             => 0,
            pbs                   => {
              "nodes"    => "1:ppn=2",
              "walltime" => "12",
              "mem"      => "26gb"
            },
          };

          my $indel_filter_level = getValue($def, "indel_filter_level", "99.7");
          my $snp_filter_level = getValue($def, "snp_filter_level", "99.7");
          my $ApplyRecalibration = "ApplyRecalibration";
          $config->{$ApplyRecalibration} = {
            class                 => "CQS::ProgramWrapperOneToOne",
            perform               => 1,
            target_dir            => "${target_dir}/" . $gatk_prefix .  getNextIndex($def, $gatk_index_snv) . "_ApplyRecalibration",
            option                => "
gatk --java-options \"-Xms5000m -Xmx6500m\" \\
  ApplyVQSR \\
  -AS \\
  -O __NAME__.indel.recalibrated.vcf.gz \\
  -V __FILE__ \\
  --recal-file __parameterFile1__ \\
  --use-allele-specific-annotations \\
  --tranches-file __parameterFile2__ \\
  --truth-sensitivity-filter-level $indel_filter_level \\
  --create-output-variant-index true \\
  -mode INDEL

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.indel.recalibrated.vcf.gz.succeed
  rm -f __NAME__.indel.recalibrated.vcf.gz.failed

  gatk --java-options \"-Xms5000m -Xmx6500m\" \\
    ApplyVQSR \\
    -O __NAME__.tmp.filtered.vcf.gz \\
    -V __NAME__.indel.recalibrated.vcf.gz \\
    --recal-file __parameterFile3__ \\
    --use-allele-specific-annotations \\
    --tranches-file __parameterFile4__ \\
    --truth-sensitivity-filter-level $snp_filter_level \\
    --create-output-variant-index true \\
    -mode SNP

  status=\$?
  if [[ \$status -eq 0 ]]; then
    touch __NAME__.filtered.vcf.gz.succeed
    rm -f __NAME__.filtered.vcf.gz.failed __NAME__.indel.recalibrated.vcf.gz __NAME__.indel.recalibrated.vcf.gz.tbi __NAME__.indel.recalibrated.vcf.gz.succeed
    mv __NAME__.tmp.filtered.vcf.gz __NAME__.filtered.vcf.gz
    mv __NAME__.tmp.filtered.vcf.gz.tbi __NAME__.filtered.vcf.gz.tbi
  else
    touch __NAME__.filtered.vcf.gz.failed
    rm -f __NAME__.filtered.vcf.gz.succeed __NAME__.tmp.filtered.vcf.gz __NAME__.tmp.filtered.vcf.gz.tbi
  fi
else
  touch __NAME__.indel.recalibrated.vcf.gz.failed
  rm -f __NAME__.indel.recalibrated.vcf.gz.succeed __NAME__.indel.recalibrated.vcf.gz __NAME__.indel.recalibrated.vcf.gz.tbi
fi
",
            interpretor           => "",
            program               => "",
            docker_prefix         => "gatk4_",
            check_program         => 0,
            source_arg            => "-V",
            source_ref            => [ $HardFilterAndMakeSitesOnlyVcf, ".variant_filtered.vcf.gz" ],
            parameterFile1_ref    => [ $IndelsVariantRecalibrator, ".indels.recal.vcf.gz" ],           
            parameterFile2_ref    => [ $IndelsVariantRecalibrator, ".indels.tranches" ],
            parameterFile3_ref    => [ $SNPsVariantRecalibrator, ".snps.recal.vcf.gz" ],           
            parameterFile4_ref    => [ $SNPsVariantRecalibrator, ".snps.tranches" ],
            output_arg            => "-O",
            output_file_prefix    => "",
            output_file_ext            => ".filtered.vcf.gz",
            output_to_same_folder => 1,
            other_localization_ext_array => [".tbi"],
            sh_direct             => 0,
            pbs                   => {
              "nodes"    => "1:ppn=1",
              "walltime" => "24",
              "mem"      => "7gb"
            },
          };

          my $FinalGatherVcf = "FinalGatherVcf";
          $config->{$FinalGatherVcf} = {
            class                 => "CQS::ProgramWrapper",
            perform               => 1,
            target_dir            => "${target_dir}/" . $gatk_prefix .  getNextIndex($def, $gatk_index_snv) . "_FinalGatherVcf",
            option                => "
set -euo pipefail

# --ignore-safety-checks makes a big performance difference so we include it in our invocation.
# This argument disables expensive checks that the file headers contain the same set of
# genotyped samples and that files are in order by position of first record.
gatk --java-options \"-Xms10g -Xmx10g\" \\
  GatherVcfsCloud \\
  --ignore-safety-checks \\
  --gather-type BLOCK \\
  --input __FILE__ \\
  --output __NAME__.tmp.vcf.gz

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.vcf.gz.succeed
  rm -f __NAME__.vcf.gz.failed 
  mv __NAME__.tmp.vcf.gz __NAME__.vcf.gz
  tabix __NAME__.vcf.gz
else
  touch __NAME__.vcf.gz.failed
  rm -f __NAME__.vcf.gz.succeed __NAME__.tmp.vcf.gz
fi
",
            interpretor           => "",
            program               => "",
            docker_prefix         => "gatk4_",
            check_program         => 0,
            source_arg            => "--input",
            source_type           => "array",
            source_ref            => [ $ApplyRecalibration, ".filtered.vcf.gz" ],
            source_join_delimiter => " \\\n  --input ",
            output_arg            => "--output",
            output_file_prefix    => "",
            output_file_ext       => ".vcf.gz",
            sh_direct             => 0,
            pbs                   => {
              "nodes"    => "1:ppn=1",
              "walltime" => "12",
              "mem"      => "10gb"
            },
          };

          my $ref_fasta_dict = getValue($def, "ref_fasta_dict");
          my $eval_interval_list = getValue($def, "eval_interval_list");
          my $CollectMetricsOnFullVcf = "CollectMetricsOnFullVcf";
          $config->{$CollectMetricsOnFullVcf} = {
            class                 => "CQS::ProgramWrapperOneToOne",
            perform               => 1,
            target_dir            => "${target_dir}/" . $gatk_prefix .  getNextIndex($def, $gatk_index_snv) . "_CollectMetricsOnFullVcf",
            option                => "
set -euo pipefail

gatk --java-options \"-Xms6000m -Xmx7000m\" \\
  CollectVariantCallingMetrics \\
  --INPUT __FILE__ \\
  --DBSNP ${dbsnp_resource_vcf} \\
  --SEQUENCE_DICTIONARY ${ref_fasta_dict} \\
  --OUTPUT __NAME__ \\
  --THREAD_COUNT 8 \\
  --TARGET_INTERVALS ${eval_interval_list}

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.succeed
  rm -f __NAME__.failed 
else
  touch __NAME__.failed
  rm -f __NAME__.succeed __NAME__.variant_calling_detail_metrics __NAME__.variant_calling_summary_metrics
fi
",
            interpretor           => "",
            program               => "",
            docker_prefix         => "gatk4_",
            check_program         => 0,
            source_arg            => "--INPUT",
            source_ref            => [ $FinalGatherVcf, ".vcf.gz" ],
            output_arg            => "--OUTPUT",
            output_file_prefix    => "",
            output_file_ext       => ".variant_calling_detail_metrics;.variant_calling_summary_metrics",
            sh_direct             => 0,
            output_to_same_folder => 1,
            pbs                   => {
              "nodes"    => "1:ppn=8",
              "walltime" => "12",
              "mem"      => "10gb"
            },
          };
          push @$summary, ($HardFilterAndMakeSitesOnlyVcf, $SitesOnlyGatherVcf, $IndelsVariantRecalibrator, $SNPsVariantRecalibrator, $ApplyRecalibration, $FinalGatherVcf, $CollectMetricsOnFullVcf);

          $merged_vcf_section = $FinalGatherVcf;
        } else {
          $merged_vcf_section = add_hard_filter_and_merge($config, $def, $summary, $target_dir, $gatk_prefix, $gatk_index_snv, $genotypeGVCFs_section);
        }
      }
    }
  }

  if(defined $merged_vcf_section){
    my $annovar_name = addAnnovar($config, $def, $summary, $target_dir, $merged_vcf_section, '.gz$', $gatk_prefix, $def, $gatk_index_snv );

    if ( $def->{annovar_param} =~ /exac/ || $def->{annovar_param} =~ /1000g/ || $def->{annovar_param} =~ /gnomad/ ) {
      my $annovar_filter_name = addAnnovarFilter( $config, $def, $summary, $target_dir, $annovar_name, $gatk_prefix, $def, $gatk_index_snv);

      if ( defined $def->{annotation_genes} ) {
        my $annovar_filter_geneannotation_name = addAnnovarFilterGeneannotation( $config, $def, $summary, $target_dir, $annovar_filter_name );
      }

      my $mafreport = addAnnovarMafReport($config, $def, $summary, $target_dir, $annovar_filter_name, $gatk_prefix, $def, $gatk_index_snv);
    }
  }

  $config->{sequencetask} = {
    class => "CQS::SequenceTaskSlurmSlim",
    perform => 1,
    target_dir            => "${target_dir}/sequencetask",
    option                => "",
    source                => {
      processing => $summary,
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
