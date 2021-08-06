#!/usr/bin/perl
package Pipeline::ExomeSeq;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use Pipeline::WdlPipeline;
use Pipeline::WGS;
use Variants::VariantsUtils;
use Data::Dumper;
use List::Util qw(max);

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performExomeSeq performExomeSeqTask addMutect2)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;
  initDefaultValue( $def, "max_thread", 8 );
  initDefaultValue( $def, "subdir",     0 );

  initDefaultValue( $def, "sra_to_fastq", 0 );

  initDefaultValue( $def, "aligner", "bwa" );
  if ( $def->{aligner} eq "bwa" ) {
    initDefaultValue( $def, "bwa_option", "" );
  }

  initDefaultValue( $def, "perform_extract_bam", 0 );
  initDefaultValue( $def, "perform_featureCounts",  0 );

  initDefaultValue( $def, "perform_gatk_callvariants",   0 );
  initDefaultValue( $def, "perform_gatk4_callvariants",  1 );
  initDefaultValue( $def, "gatk4_scatter_count",  0 );
  initDefaultValue( $def, "callvariants_vqsr_mode", 1 );
  initDefaultValue( $def, "gatk_callvariants_vqsr_mode", getValue($def, "callvariants_vqsr_mode") );
  initDefaultValue( $def, "has_chr_in_chromosome_name" , 0);
  initDefaultValue( $def, "perform_gatk4_pairedfastq2bam",  0 );

  initDefaultValue( $def, "perform_target_coverage",  0 );

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
  initDefaultValue( $def, "perform_muTect2",       0 );
  initDefaultValue( $def, "perform_annovar",      0 );
  initDefaultValue( $def, "perform_cnv",          1 );
  initDefaultValue( $def, "perform_vep",          0 );
  initDefaultValue( $def, "perform_IBS",          0 );

  if ( $def->{perform_muTect} || $def->{perform_muTect2indel} || $def->{perform_muTect2} ) {
    initDefaultValue( $def, "perform_CrosscheckFingerprints", defined $def->{hapmap_file} );
    if ( defined $def->{mills} ) {
      initDefaultValue( $def, "indel_realignment", 1 );
      initDefaultValue( $def, "indel_vcf_files",   $def->{mills} );
    }
    else {
      initDefaultValue( $def, "indel_realignment", 0 );
    }
  }else{
    initDefaultValue( $def, "perform_CrosscheckFingerprints", 0 );
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

sub addMutect2 {
  my ($config, $def, $tasks, $target_dir, $bam_input, $mutect2_call, $option, $use_germline_resource, $pon) = @_;

  my $run_funcotator;
  if(defined $def->{run_funcotator}){
    $run_funcotator = $def->{run_funcotator};
  }else{
    if ($def->{ncbi_build} eq "GRCh38") { #based on genome, hg38=true, else false
      $run_funcotator="true";
    }else{
      $run_funcotator="false";
    }
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

  my $fasta_file = $def->{ref_fasta};
  if(not defined $fasta_file){
    $fasta_file = $def->{fasta_file};
  }
  die "define ref_fasta or fasta_file for mutect2" if not defined $fasta_file;

  $config->{$mutect2_call} = {     
    "class" => "GATK4::MuTect2",
    "target_dir" => "${target_dir}/$mutect2_call",
    "option" => $option,
    "m2_extra_filtering_args" => $def->{m2_extra_filtering_args},
    "source_ref" => [ $bam_input, '.bam$' ],
    "variants_for_contamination" => $def->{variants_for_contamination},
    "run_orientation_bias_mixture_model_filter" => $def->{'Mutect2.run_orientation_bias_mixture_model_filter'},
    "fasta_file" => $fasta_file,
    pbs=> {
      "nodes"     => "1:ppn=1",
      "walltime"  => getValue($def, "mutect2_walltime", "22"),
      "mem"       => getValue($def, "mutect2_memory", "40gb"),
    },
  };

  if(defined $def->{covered_bed}) {
    $config->{$mutect2_call}{intervals} = $def->{covered_bed};
    #print($def->{covered_bed});
  }elsif(defined $config->{intervals_ref}){
    $config->{$mutect2_call}{intervals_ref} = $config->{intervals_ref};
  }

  if($use_germline_resource){
    $config->{$mutect2_call}{germline_resource} = $def->{germline_resource};
  }

  if(defined $pon){
    my $pon_key = (-e $pon) ? "panel_of_normals" : "panel_of_normals_ref";
    $config->{$mutect2_call}{$pon_key} = $pon;
  }

  push @$tasks, $mutect2_call;
}

sub init_muTect2_groups {
  my ($config, $def) = @_;

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
}

sub add_muTect2_PON {
  my ($config, $def, $individual, $summary, $target_dir, $bam_input, $task_prefix) = @_;
    
  my $mutect2_pon_normal_files="mutect2_pon_normal_files";

  if($def->{pon_samples}){
    #https://gatk.broadinstitute.org/hc/en-us/articles/360035531132
    $config->{$mutect2_pon_normal_files} = {     
      "class" => "CQS::SamplePickTask",
      "source_ref" => $bam_input,
      "sample_names" => $def->{pon_samples},
    };
  }else{
    init_muTect2_groups($config, $def);

    if(is_muTect2_tumor_only($config)){
      die "Cannot generate PON for tumor only project!";
    }

    $config->{$mutect2_pon_normal_files} = {     
      "class" => "CQS::GroupPickTask",
      "source_ref" => $bam_input,
      "groups_ref" => "mutect2_groups",
      "sample_index_in_group" => 0, 
    };
  }
  
  my $perform_mutect2_by_wdl = getValue($def, "perform_mutect2_by_wdl", 1);

  my $pon_mutect2call_ref;
  my $pon_mutect2_option = getValue($def, "pon_mutect2_option", "--max-mnp-distance 0");
  if($perform_mutect2_by_wdl){
    my $pon_mutect2call = addMutect2Wdl($config, $def, $individual, $target_dir, $task_prefix, $pon_mutect2_option, 1, 1, $mutect2_pon_normal_files, undef, undef);
    $pon_mutect2call_ref = [ $pon_mutect2call, '.vcf$' ];
  }else{
    my $pon_mutect2call = $task_prefix . "01_call";
    addMutect2($config, $def, $individual, $target_dir, $mutect2_pon_normal_files, $pon_mutect2call, $pon_mutect2_option, 0, undef);
    $pon_mutect2call_ref = [ $pon_mutect2call, '.pass.vcf.gz$' ];
  }

  my $pon_task = $task_prefix . "02_table";
  $config->{$pon_task} = {
    class => "GATK4::CreateSomaticPanelOfNormals",
    option => "",
    source_ref => $pon_mutect2call_ref,
    target_dir   => "${target_dir}/$pon_task",
    germline_resource => $def->{germline_resource},
    "intervals" => $def->{covered_bed},
    "fasta_file" => $def->{ref_fasta},
    sh_direct    => 1,
    pbs          => {
      "nodes"    => "1:ppn=1",
      "walltime" => getValue($def, "pon_walltime", "22"),
      "mem"      => getValue($def, "pon_memory", "40gb"),
    },
  };
  push (@$summary, $pon_task);
  return($pon_task);
}

sub init_normal_tumor_files {
  my ($config, $def, $bam_input) = @_;

  if($def->{pon_samples}){
    #https://gatk.broadinstitute.org/hc/en-us/articles/360035531132
    $config->{"muTect2_tumor_files"} = {     
      "class" => "CQS::SamplePickTask",
      "source_ref" => $bam_input,
      "not_sample_names" => $def->{pon_samples},
    };
    return("muTect2_tumor_files");
  }else{
    init_muTect2_groups($config, $def);

    my $is_tumor_only = is_muTect2_tumor_only($config);
    if($is_tumor_only){
      $config->{"muTect2_tumor_files"} = {     
        "class" => "CQS::GroupPickTask",
        "source_ref" => $bam_input,
        "groups_ref" => "mutect2_groups",
        "sample_index_in_group" => 0, 
      };
      return("muTect2_tumor_files");
    }else{
      $config->{"muTect2_tumor_files"} = {     
        "class" => "CQS::GroupPickTask",
        "source_ref" => $bam_input,
        "groups_ref" => "mutect2_groups",
        "sample_index_in_group" => 1, 
      };
      $config->{"muTect2_normal_files"} = {     
        "class" => "CQS::GroupPickTask",
        "source_ref" => $bam_input,
        "groups_ref" => "mutect2_groups",
        "sample_index_in_group" => 0, 
      };
      return("muTect2_tumor_files", "muTect2_normal_files");
    }
  }
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  if($def->{files} && $def->{groups}){
    my $groups = $def->{groups};
    my $files = $def->{files};
    my @errors;
    for my $gname (keys %$groups){
      my $samples = $groups->{$gname};
      for my $sample (@$samples){
        if (not defined $files->{$sample}){
          push(@errors, "ERROR: $sample in group $gname was not defined in files.");
        }
      }
    }
    if(scalar(@errors) > 0){
      for my $error (@errors){
        print $error . "\n";
      }
      die "Please have a check at configuration.";
    }
  }

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  if ($def->{perform_preprocessing_only}) {
    return $config;
  }

  my $step3 = [];
  my $step4 = [];
  my $step5 = [];
  my $step6 = [];

  my $email      = getValue( $def, "email" );
  my $max_thread = getValue( $def, "max_thread" );
  my $covered_bed = getValue( $def, "covered_bed" );

  my $geneLocus = undef;
  my $chrCode = getValue($def, "has_chr_in_chromosome_name") ? ";addChr=1" : "";

  my $species = getValue($def, "species");

  if ( defined $def->{annotation_genes} ) {
    $geneLocus = "annotation_genes_locus";
    $config->{$geneLocus} = {
      class      => "CQS::UniqueR",
      perform    => 1,
      target_dir => $target_dir . '/' . $geneLocus,
      rtemplate  => "../Annotation/getGeneLocus.r",
      rCode      => "host=\""
        . getValue( $def, "biomart_host" )
        . "\";dataset=\""
        . getValue( $def, "biomart_dataset" )
        . "\";symbolKey=\""
        . getValue( $def, "biomart_symbolKey" )
        . "\";genesStr=\""
        . getValue( $def, "annotation_genes" ) . "\"" . $chrCode,
      output_file_ext => ".missing;.bed",
      sh_direct       => 1,
      'pbs'           => {
        'nodes'    => '1:ppn=1',
        'mem'      => '40gb',
        'walltime' => '10'
      },
    };
    push( @$summary, $geneLocus );
  }

  my $bam_ref;
  my $bam_input;
  my $fasta;
  my $dbsnp = $def->{dbsnp};
  my $mills = $def->{mills};
  my $gatk_jar   = getValue( $def, "gatk3_jar" );
  my $picard_jar = getValue( $def, "picard_jar" );

  my $rg_name_regex = undef;
  my $alignment_source_ref = $source_ref;
  if ($def->{perform_gatk4_pairedfastq2bam}){
    $bam_input = addPairedFastqToProcessedBam($config, $def, $individual, $target_dir, $alignment_source_ref);
    $bam_ref = [$bam_input, ".bam\$"];
    $fasta = getValue( $def, "bwa_fasta" );

    #TEQC for target region coverage
    $config->{"TEQC"} = {
      class      => "CQS::UniqueR",
      perform    => 1,
      target_dir => $target_dir . '/' . "TEQC",
      parameterSampleFile1_ref=> $bam_ref,
      parameterFile1=> $covered_bed,
      rtemplate  => "runTEQC.R",
      rCode      => "genome=\""
        . getValue($def, "annovar_buildver", "hg38")
        . "\";",
      output_file_ext => ".TEQC",
      sh_direct       => 1,
      'pbs'           => {
        'nodes'    => '1:ppn=1',
        'mem'      => '20gb',
        'walltime' => '10'
      },
    };

  }else{
    if ($def->{aligner_scatter_count}){
      my $splitFastq = "bwa_01_splitFastq";
      $config->{ $splitFastq } = {
        class       => "CQS::ProgramWrapperOneToMany",
        perform     => 1,
        target_dir  => "$target_dir/$splitFastq",
        option      => "",
        interpretor => "python3",
        program     => "../Format/splitFastq.py",
        source_ref      => $source_ref,
        source_arg            => "-i",
        source_join_delimiter => ",",
        output_to_same_folder => 1,
        output_arg            => "-o",
        output_file_prefix    => "",
        output_file_ext       => "._ITER_.1.fastq.gz",
        output_other_ext      => "._ITER_.2.fastq.gz",
        iteration_arg         => "--trunk",
        iteration             => $def->{aligner_scatter_count},
        sh_direct             => 1,
        pbs                   => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "10",
          "mem"       => "10gb"
        },
      };
      $rg_name_regex = "(.+)_ITER_";
      $alignment_source_ref = $splitFastq;
      push @$individual, $splitFastq;
    }

    #based on paper https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1097-3, we don't do markduplicate anymore
    if ( $def->{aligner} eq "bwa") {
      $fasta = getValue( $def, "bwa_fasta" );
      my $bwa = $def->{aligner_scatter_count}?"bwa_02_alignment":"bwa";
      $config->{ $bwa } = {
        class                 => "Alignment::BWA",
        perform               => 1,
        target_dir            => "${target_dir}/" . getNextFolderIndex($def) . $bwa,
        option                => getValue( $def, "bwa_option" ),
        bwa_index             => $fasta,
        source_ref            => $alignment_source_ref,
        use_tmp_folder        => $def->{"bwa_use_tmp_folder"},
        use_sambamba          => $def->{"use_sambamba"},
        output_to_same_folder => 1,
        rg_name_regex         => $rg_name_regex,
        sh_direct             => 0,
        pbs                   => {
          "nodes"    => "1:ppn=" . $max_thread,
          "walltime" => getValue($def, "bwa_walltime", "22"),
          "mem"      => getValue($def, "bwa_memory", "40gb")
        },
      };
      $bam_ref = [ $bwa, ".bam\$" ];
      $bam_input = $bwa;
      push @$individual, ( $bwa );

      my $bwa_summary = $def->{aligner_scatter_count}?"bwa_04_summary":"bwa_summary";
      if ($def->{aligner_scatter_count}) {
        add_BWA_summary($config, $def, $summary, $target_dir, $bwa_summary, $bwa, $rg_name_regex);
      }else{
        add_BWA_summary($config, $def, $summary, $target_dir, $bwa_summary, $bwa);
      }
    }
    else {
      die "Unknown alinger " . $def->{aligner};
    }

    if ($def->{aligner_scatter_count}){
      my $mergeBam = "bwa_03_merge";
      $config->{$mergeBam} = {
        class                 => "CQS::ProgramWrapperManyToOne",
        perform               => 1,
        target_dir            => "$target_dir/$mergeBam",
        option                => "-t 8 __OUTPUT__ __INPUT__",
        interpretor           => "",
        check_program         => 0,
        program               => "sambamba merge",
        source_ref            => $bam_ref,
        source_arg            => "",
        source_join_delimiter => " ",
        output_to_same_folder => 1,
        output_arg            => "-o",
        output_file_prefix    => ".bam",
        output_file_ext       => ".bam",
        sh_direct             => 1,
        pbs                   => {
          "nodes"     => "1:ppn=8",
          "walltime"  => "10",
          "mem"       => "10gb"
        },
      };
      $bam_input = $mergeBam;
      $bam_ref = [ $mergeBam, ".bam\$" ];
      push @$individual, (  $mergeBam );
    }

    if(getValue($def, "perform_bam_validation", 0)){
      add_bam_validation($config, $def, $individual, $target_dir, $bam_input . "_bam_validation", $bam_ref );
    }

    my $perform_cnv = $def->{perform_cnv_cnMOPs} || $def->{perform_cnv_gatk4_cohort} || $def->{perform_cnv_xhmm};

    my $indel_vcf_files;
    if ( $def->{indel_realignment} ) {
      if ( !defined $def->{indel_vcf_files} & defined $dbsnp ) {
        $indel_vcf_files = $dbsnp;
      }
      else {
        $indel_vcf_files = getValue( $def, "indel_vcf_files" );
      }
    }
    my $vcf;
    if ( defined $dbsnp & defined $mills ) {
      $vcf = [ $dbsnp, $mills ];
    }
    elsif ( defined $dbsnp ) {
      $vcf = [$dbsnp];
    }
    elsif ( defined $mills ) {
      $vcf = [$mills];
    }
    else {
      $vcf = undef;
    }

    #based on paper https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1097-3, we don't do markduplicate anymore
    my $refine_name = getValue($def, "perform_gatk4_refine", 0) ? $def->{aligner} . "_g4_refine" : $def->{aligner} . "_refine";
    my $refine_memory = getValue($def, "refine_memory", "40gb");
    $config->{$refine_name} = {
      class      => getValue($def, "perform_gatk4_refine", 0) ? "GATK4::Refine" : "GATK::Refine",
      perform    => 1,
      target_dir => "${target_dir}/$refine_name",
      option     => "-Xmx40g",

      #gatk_option => "--fix_misencoded_quality_scores",
      gatk_option              => "",
      fasta_file               => $fasta,
      source_ref               => $bam_ref,
      vcf_files                => $vcf,
      gatk_jar                 => $gatk_jar,
      picard_jar               => $picard_jar,
      remove_duplicate         => getValue($def, "remove_duplicate"),
      use_tmp_folder           => $def->{"bwa_use_tmp_folder"},
      sh_direct                => 0,
      slim_print_reads         => 1,
      samtools_baq_calibration => 0,
      indel_realignment        => $def->{indel_realignment},
      indel_vcf_files          => $indel_vcf_files,
      sorted                   => 1,
      pbs                      => {
        "nodes"    => "1:ppn=1",
        "walltime" => "48",
        "mem"      => $refine_memory
      },
    };
    push @$individual, ($refine_name);

    $bam_input = $refine_name;
    $bam_ref = [ $refine_name, ".bam\$" ];

    add_alignment_summary($config, $def, $summary, $target_dir, "${refine_name}_summary", "../Alignment/AlignmentUtils.r;../Alignment/BWASummary.r", ".chromosome.csv;.chromosome.png", undef, [$refine_name, ".chromosome.count"] );

    if(getValue($def, "perform_bam_validation", 0)){
      add_bam_validation($config, $def, $individual, $target_dir, $refine_name . "_bam_validation", $bam_ref );
    }
  }

  if($def->{filter_soft_clip} && ((!defined $def->{perform_gatk4_pairedfastq2bam}) || (!$def->{perform_gatk4_pairedfastq2bam}))){
    my $soft_clip_name = $bam_input . "_nosoftclip";
    $config->{$soft_clip_name} = {
      class                 => "CQS::ProgramWrapperOneToOne",
      perform               => 1,
      target_dir            => "${target_dir}/${soft_clip_name}",
      option                => "--min-mapq " . getValue($def, "soft_clip_min_mapq", 10) . " -i __FILE__ -o __NAME__.nosoftclip.bam

samtools index __NAME__.nosoftclip.bam.bai
samtools idxstats __NAME__.nosoftclip.bam > __NAME__.nosoftclip.bam.chromosome.count
",
      interpretor           => "python3",
      program               => "../GATK/filterSoftClip.py",
      source_arg            => "-i",
      source_ref            => $bam_ref,
      output_to_same_folder => 1,
      use_tmp_folder        => 1,
      output_arg            => "-o",
      output_file_prefix    => ".nosoftclip.bam",
      output_file_ext       => ".nosoftclip.bam;.nosoftclip.bam.bai;.nosoftclip.bam.chromosome.count",
      use_tmp_folder        => 1,
      sh_direct             => 0,
      pbs                   => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "10gb"
      },
    };
    push @$individual, ($soft_clip_name);
    $bam_input = $soft_clip_name;
    $bam_ref = [$bam_input, ".bam\$"];

    add_alignment_summary($config, $def, $summary, $target_dir, "${soft_clip_name}_summary", "../Alignment/AlignmentUtils.r;../Alignment/BWASummary.r", ".chromosome.csv;.chromosome.png", undef, [$soft_clip_name, ".chromosome.count"] );
  }

  if ($def->{perform_TEQC}) {
    $fasta = getValue( $def, "bwa_fasta" );

    #TEQC for target region coverage
    $config->{"TEQC"} = {
      class      => "CQS::UniqueR",
      perform    => 1,
      target_dir => $target_dir . '/' . "TEQC",
      parameterSampleFile1_ref=> $bam_ref,
      parameterFile1=> $covered_bed,
      rtemplate  => "runTEQC.R",
      rCode      => "genome=\""
        . getValue($def, "annovar_buildver", "hg38")
        . "\";",
      output_file_ext => ".TEQC",
      sh_direct       => 1,
      'pbs'           => {
        'nodes'    => '1:ppn=1',
        'mem'      => '20gb',
        'walltime' => '10'
      },
    };
    push @$summary, ("TEQC");
  }

  if ($def->{perform_CrosscheckFingerprints}){
      my $CrosscheckFingerprints_name = $bam_input . "_CrosscheckFingerprints";
      my $files = $def->{files};
      my $files_count = scalar(keys %$files);
      my $height = max($files_count * 60, 3000);
      my $width = $height + 300;
      my $hapmap = getValue($def, "hapmap_file");
      $config->{$CrosscheckFingerprints_name} = {
        class                 => "CQS::ProgramWrapper",
        perform               => 1,
        target_dir            => "${target_dir}/${CrosscheckFingerprints_name}",
        option                => "
if [[ -e __NAME__.failed ]]; then
  rm __NAME__.failed
fi

java -Xmx40g -jar $picard_jar CrosscheckFingerprints INPUT=__FILE__ \\
  H=$hapmap \\
  CROSSCHECK_BY=SAMPLE \\
  LOD_THRESHOLD=-5 \\
  EXPECT_ALL_GROUPS_TO_MATCH=true \\
  NUM_THREADS=8 \\
  OUTPUT=__NAME__.crosscheck_metrics 

status=\$?
if [[ \$status -ne 0 ]]; then
  touch __NAME__.failed
else
  java -Xmx40g -jar $picard_jar ClusterCrosscheckMetrics INPUT=__NAME__.crosscheck_metrics \\
    LOD_THRESHOLD=-5 \\
    OUTPUT=__NAME__.clustered.crosscheck_metrics

  cat >__NAME__.r <<EOT
library('ggplot2')
mat<-read.table('__NAME__.crosscheck_metrics', sep='\\t', header=T)

png('__NAME__.crosscheck_metrics.png', width=$width, height=$height, res=300)

g<-ggplot(mat, aes(x=LEFT_GROUP_VALUE, y=RIGHT_GROUP_VALUE)) + 
  geom_point(aes(size=LOD_SCORE, color=RESULT)) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0),
        axis.title = element_blank())
print(g)
dev.off()

EOT

  R --vanilla -f __NAME__.clustered.crosscheck_metrics.r
fi

",
        interpretor           => "",
        check_program         => 0,
        program               => "",
        parameterSampleFile1_arg    => "INPUT=",
        parameterSampleFile1_ref    => $bam_ref,
        parameterSampleFile1_fileonly    => 1,
        output_to_same_folder => 1,
        no_output             => 1,
        output_arg            => "OUTPUT=",
        output_file_ext       => ".clustered.crosscheck_metrics",
        sh_direct             => 1,
        pbs                   => {
          "nodes"     => "1:ppn=8",
          "walltime"  => "10",
          "mem"       => "40gb"
        },
      };
      push @$summary, $CrosscheckFingerprints_name;
  }

  if ($def->{perform_target_coverage}){
    my $target_coverage_task = $bam_input . "_target_coverage";
    my $bait_intervals = getValue($def, "bait_intervals_file");
    my $target_intervals = getValue($def, "target_intervals_file");
    $config->{$target_coverage_task} = {
      class                 => "CQS::ProgramWrapperOneToOne",
      perform               => 1,
      target_dir            => "${target_dir}/${target_coverage_task}",
      option                => "--java-options \"-Xms__MEMORY__\" CollectHsMetrics \\
      --INPUT __FILE__ \\
      --OUTPUT __OUTPUT__ \\
      --BAIT_INTERVALS ${bait_intervals} \\
      --TARGET_INTERVALS ${target_intervals}
  ",
      interpretor           => "",
      docker_prefix         => "gatk4_",
      program               => "gatk",
      check_program         => 0,
      source_arg            => "",
      source_ref            => $bam_ref,
      source_join_delimiter => " ",
      output_to_same_folder => 1,
      output_arg            => "",
      output_file_prefix    => "_hs_metrics.txt",
      output_file_ext       => "_hs_metrics.txt",
      sh_direct             => 0,
      pbs                   => {
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    };
    push @$individual, ($target_coverage_task);
  }

  if($def->{perform_bamsnap} && $def->{"bamsnap_locus"}){
    #addBamsnapLocus($config, $def, $summary, $target_dir, "bamsnap_locus", $bam_ref);
    addBamsnapLocus($config, $def, $summary, $target_dir, "bamsnap_locus", ['bwa', '.bam$']);
  }

  if($def->{perform_extract_bam}){
    my $extract_bam_locus = getValue($def, "extract_bam_locus");
    my $extract_bam_task = "extract_bam_locus";
    $config->{$extract_bam_task} = {
      class => "CQS::ProgramWrapperOneToOne",
      target_dir => $target_dir . "/" . getNextFolderIndex($def) . $extract_bam_task,
      interpretor => "",
      check_program => 0,
      option => "view -b -o __OUTPUT__ __FILE__ " . $extract_bam_locus . "; samtools index __OUTPUT__; ",
      program => "samtools",
      source_ref => [ $bam_input, ".bam\$" ],
      output_arg => "",
      output_file_prefix => ".bam",
      output_file_ext => ".bam",
      output_to_same_folder => 1,
      sh_direct   => 1,
      pbs => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "10gb"
      }
    };
    push @$individual, ($extract_bam_task);
  }

  if($def->{perform_featureCounts}){
    my $featureCounts = $bam_input . "_featureCounts";
    my $featureCountFolder = "${target_dir}/${featureCounts}";
    $config->{$featureCounts} = {
      class      => "Count::FeatureCounts",
      perform    => 1,
      target_dir => $featureCountFolder,
      option     => "-F SAF",
      source_ref => $bam_input,
      gff_file   => getValue($def, "saf_file"),
      is_paired_end  => is_paired_end($def),
      sh_direct  => 1,
      pbs        => {
        "email"     => $email,
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "23",
        "mem"       => "40gb"
      },
    };
    push @$individual, ($featureCounts);

    my $featureCountsSummary = $featureCounts . "_summary";
    $config->{$featureCountsSummary} = {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      target_dir               => $featureCountFolder,
      option                   => "",
      rtemplate                => "../Alignment/STARFeatureCount.r",
      output_file_ext          => ".FeatureCountSummary.csv",
      output_other_ext         => ".FeatureCountSummary.csv.png",
      parameterSampleFile2_ref => [ $featureCounts, ".count.summary" ],
      sh_direct                => 1,
      pbs                      => {
        "email"     => $email,
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "2",
        "mem"       => "10gb"
      },
    };

    push @$summary, $featureCountsSummary;

    my $name_map_file = $def->{name_map_file};
    my $countTable = $featureCounts . "_table";
    $config->{$countTable} = {
      class                     => "CQS::CQSDatatable",
      perform                   => 1,
      target_dir                => $target_dir . "/" . $countTable,
      option                    => "-k 0 -v 6 -e --fillMissingWithZero",
      source_ref                => [$featureCounts, ".count\$"],
      output_proteincoding_gene => 0,
      name_map_file             => $name_map_file,
      sh_direct                 => 1,
      pbs                       => {
        "email"     => $email,
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "23",
        "mem"       => "10gb"
      },
    };

    push @$summary, $countTable;
  }

  my $gatk_index = $def;
  my $gatk_index_snv = "SNV_Index";

  my $gatk_prefix = "";
  my $snv_index = 0;
  my $filter_name = "";
  if ( $def->{perform_gatk4_callvariants} ) {
    $gatk_prefix = $bam_input . "_gatk4_SNV_";

    my $perform_gatk4_by_scatter = getValue($def, "perform_gatk4_by_scatter", 0);

    if ($perform_gatk4_by_scatter){
      if (not defined $def->{interval_list_file}){
        my $gatk4_scatter_count = getValue($def, "gatk4_scatter_count", 0);
        my $scatterIntervalTask = $gatk_prefix . "00_scatter_intervals";
        $config->{$scatterIntervalTask} = {
          class                 => "CQS::ProgramWrapperOneToOne",
          perform               => 1,
          target_dir            => "${target_dir}${scatterIntervalTask}",
          option                => "--java-options \"-Xmx40g\" SplitIntervals \\
          -L __FILE__ \\
          -O __NAME__.intervals \\
          -scatter $gatk4_scatter_count \\
          -R $fasta 

ls \$(pwd)/__NAME__.intervals/* > __NAME__.intervals_list
",
          interpretor           => "",
          docker_prefix         => "gatk4_",
          program               => "gatk",
          check_program         => 0,
          source_arg            => "",
          source            => {
            $def->{task_name} => [$covered_bed]
          },
          output_to_same_folder => 1,
          output_arg            => "",
          no_output => 1,
          output_file_prefix    => ".intervals_list",
          output_file_ext       => ".intervals_list",
          sh_direct             => 1,
          pbs                   => {
            "nodes"    => "1:ppn=1",
            "walltime" => "2",
            "mem"      => "10gb"
          },
        };

        performTask($config, $scatterIntervalTask);
        die "\n\nRun task ${target_dir}${scatterIntervalTask} first, then add following line to configuration: \n\ninterval_list_file => '" . $config->{$scatterIntervalTask}{target_dir} . "/" . $def->{task_name} . ".intervals_list', \n\n";
      }
    }

    my $gvcf_name         = $gatk_prefix . getNextIndex($gatk_index, $gatk_index_snv) . "_hc_gvcf";
    $config->{$gvcf_name} = {
      class             => "GATK4::HaplotypeCaller",
      perform           => 1,
      target_dir        => "${target_dir}/$gvcf_name",
      option            => "",
      source_ref        => $bam_ref,
      java_option       => "",
      fasta_file        => $fasta,
      extension         => ".g.vcf.gz",
      bed_file          => $covered_bed,
      blacklist_file    => $def->{blacklist_file},
      by_chromosome     => 0,
      gvcf              => 1,
      sh_direct         => 0,
      pbs               => {
        "email"    => $email,
        "nodes"    => "1:ppn=" . $max_thread,
        "walltime" => "24",
        "mem"      => "40gb"
      },
    };
    push @$individual, ($gvcf_name);

    if(not getValue($def, "callvariants_vqsr_mode")){
      my $genotypeGVCFs_section;
      
      if ($perform_gatk4_by_scatter) {
        $genotypeGVCFs_section = add_gvcf_to_genotype_scatter($config, $def, $summary, $target_dir, $gatk_prefix, $gatk_index_snv, $gvcf_name, "covered_bed");
      }else{
        $genotypeGVCFs_section = add_gvcf_to_genotype($config, $def, $summary, $target_dir, $gatk_prefix, $gatk_index_snv, $gvcf_name, "covered_bed");
      }
      $filter_name = add_hard_filter_and_merge($config, $def, $summary, $target_dir, $gatk_prefix, $gatk_index_snv, $genotypeGVCFs_section);
    }
    elsif(getValue($def, "gatk4_variant_filter_by_chromosome", 0)){
      my $vqsr_prefix = $gatk_prefix . getNextIndex($gatk_index, $gatk_index_snv) ;
      my $filter_name_chr = $vqsr_prefix . "_vqsr_1_scatter";
      $config->{$filter_name_chr} = {
        class             => "GATK4::VariantFilterScatterChromosome",
        perform           => 1,
        target_dir        => "${target_dir}/$filter_name_chr",
        option            => "",
        vqsr_mode         => 1,
        source_ref        => "$gvcf_name",
        java_option       => "",
        fasta_file        => $fasta,
        dbsnp_vcf         => $dbsnp,
        chromosome_names  => getValue($def, "chromosome_names"),
        sh_direct         => 0,
        pbs               => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "4",
          "mem"      => "40gb"
        },
      };
      push @$summary, ($filter_name_chr);

      my $filter_name_chr_recalibrator = $vqsr_prefix . "_vqsr_2_recalibrator";
      $config->{$filter_name_chr_recalibrator} = {
        class             => "GATK4::VariantFilterRecalibrator",
        perform           => 1,
        target_dir        => "${target_dir}/$filter_name_chr_recalibrator",
        option            => "",
        vqsr_mode         => 1,
        source_ref        => ["$filter_name_chr",  "sites_only.vcf.gz"],
        java_option       => "",
        fasta_file        => $fasta,
        dbsnp_vcf         => $dbsnp,
        hapmap_vcf        => $def->{hapmap},
        omni_vcf          => $def->{omni},
        g1000_vcf         => $def->{g1000},
        axiomPoly_vcf     => $def->{axiomPoly},
        mills_vcf         => $mills,
        chromosome_names  => getValue($def, "chromosome_names"),
        sh_direct         => 0,
        pbs               => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "4",
          "mem"      => "10gb"
        },
      };
      push @$summary, ($filter_name_chr_recalibrator);

      my $filter_name_chr_recalibrator_apply = $vqsr_prefix . "_vqsr_3_applyVQSR";
      $config->{$filter_name_chr_recalibrator_apply} = {
        class             => "GATK4::VariantFilterApplyVQSR",
        perform           => 1,
        target_dir        => "${target_dir}/$filter_name_chr_recalibrator_apply",
        option            => "",
        vqsr_mode         => 1,
        source_ref        => [$filter_name_chr,  "variant_filtered.vcf.gz"],
        java_option       => "",
        indels_recalibration_ref => [$filter_name_chr_recalibrator, ".indels.recal.vcf.gz"],
        indels_tranches_ref => [$filter_name_chr_recalibrator, ".indels.tranches"],
        snps_recalibration_ref => [$filter_name_chr_recalibrator, ".snp.recal.vcf.gz"],
        snps_tranches_ref => [$filter_name_chr_recalibrator, ".snp.tranches"],
        chromosome_names  => getValue($def, "chromosome_names"),
        sh_direct         => 1,
        pbs               => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "4",
          "mem"      => "10gb"
        },
      };
      push @$summary, ($filter_name_chr_recalibrator_apply);

      my $filter_name_chr_recalibrator_apply_gather = $vqsr_prefix . "_vqsr_4_gather";
      $config->{$filter_name_chr_recalibrator_apply_gather} = {
        class             => "GATK4::VariantFilterGather",
        perform           => 1,
        target_dir        => "${target_dir}/$filter_name_chr_recalibrator_apply_gather",
        option            => "",
        source_ref        => ["$filter_name_chr_recalibrator_apply",  "pass.vcf.gz"],
        fasta_file        => $fasta,
        java_option       => "",
        chromosome_names  => getValue($def, "chromosome_names"),
        sh_direct         => 1,
        pbs               => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "4",
          "mem"      => "10gb"
        },
      };
      push @$summary, ($filter_name_chr_recalibrator_apply_gather);

      $filter_name = $filter_name_chr_recalibrator_apply_gather;
    }else{
      $filter_name = $gatk_prefix . getNextIndex($gatk_index, $gatk_index_snv) . "_vqsr";
      $config->{$filter_name} = {
        class             => "GATK4::VariantFilter",
        perform           => 1,
        target_dir        => "${target_dir}/$filter_name",
        option            => "",
        vqsr_mode         => 1,
        source_ref        => "$gvcf_name",
        java_option       => "",
        fasta_file        => $fasta,
        dbsnp_vcf         => $dbsnp,
        hapmap_vcf        => $def->{hapmap},
        omni_vcf          => $def->{omni},
        g1000_vcf         => $def->{g1000},
        axiomPoly_vcf     => $def->{axiomPoly},
        mills_vcf         => $mills,
        species    => $species,
        sh_direct         => 1,
        pbs               => {
          "email"    => $email,
          "nodes"    => "1:ppn=8",
          "walltime" => "24",
          "mem"      => "40gb"
        },
      };
      push @$summary, ($filter_name);
    }
  }
  elsif ( $def->{perform_gatk_callvariants} ) {
    $gatk_prefix = $bam_input . "_gatk3_SNV_Germline_";
    my $gvcf_name = $gatk_prefix . getNextIndex($gatk_index, $gatk_index_snv) . "_hc_gvcf";
    $config->{$gvcf_name} = {
      class         => "GATK::HaplotypeCaller",
      perform       => 1,
      target_dir    => "${target_dir}/$gvcf_name",
      option        => "",
      source_ref    => $bam_ref,
      java_option   => "",
      fasta_file    => $fasta,
      gatk_jar      => $gatk_jar,
      extension     => ".g.vcf",
      bed_file      => $covered_bed,
      by_chromosome => 0,
      gvcf          => 1,
      sh_direct     => 0,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=" . $max_thread,
        "walltime" => "24",
        "mem"      => "40gb"
      },
    };
    push @$individual, ($gvcf_name);

    if ( $def->{gatk_callvariants_vqsr_mode} ) {
      $filter_name = $gatk_prefix . getNextIndex($gatk_index, $gatk_index_snv) . "_vqsr";
      $config->{$filter_name} = {
        class       => "GATK::VariantFilter",
        perform     => 1,
        target_dir  => "${target_dir}/$filter_name",
        option      => "",
        vqsr_mode   => 1,
        source_ref  => "$gvcf_name",
        java_option => "",
        fasta_file  => $fasta,
        dbsnp_vcf   => $dbsnp,
        hapmap_vcf  => $def->{hapmap},
        omni_vcf    => $def->{omni},
        g1000_vcf   => $def->{g1000},
        mills_vcf   => $mills,
        gatk_jar    => $gatk_jar,
        sh_direct   => 1,
        pbs         => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "40gb"
        },
      };
    }
    else {
      $filter_name = $gatk_prefix . getNextIndex($gatk_index, $gatk_index_snv) . "_hardfilter";

      $config->{$filter_name} = {
        class       => "GATK::VariantFilter",
        perform     => 1,
        target_dir  => "${target_dir}/$filter_name",
        option      => "",
        source_ref  => $gvcf_name,
        java_option => "",
        gatk_jar    => $gatk_jar,
        fasta_file  => $fasta,
        sh_direct   => 1,
        vqsr_mode   => 0,
        is_rna      => 0,
        pbs         => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "40gb"
        },
      };
    }
    push @$summary, ($filter_name);
  }

  #print($filter_name);

  my $annovar_filter_geneannotation_name = undef;
  if ( $def->{perform_gatk4_callvariants} or $def->{perform_gatk_callvariants} ) {
    if ( $def->{filter_variants_by_allele_frequency} ) {
      my $maf_filter_name = $gatk_prefix . getNextIndex($gatk_index, $gatk_index_snv) . "_filterMAF";
      add_maf_filter($config, $def, $summary, $target_dir, $maf_filter_name, $filter_name);
      $filter_name = $maf_filter_name;
    }

    if ( $def->{perform_annovar} ) {
      my $annovar_name = addAnnovar( $config, $def, $summary, $target_dir, $filter_name, undef, $gatk_prefix, $gatk_index, $gatk_index_snv );

      if ( $def->{annovar_param} =~ /exac/ || $def->{annovar_param} =~ /1000g/ || $def->{annovar_param} =~ /gnomad/ ) {
        my $annovar_filter_name = addAnnovarFilter( $config, $def, $summary, $target_dir, $annovar_name, $gatk_prefix, $gatk_index, $gatk_index_snv);

        if ( defined $def->{annotation_genes} ) {
          $annovar_filter_geneannotation_name = addAnnovarFilterGeneannotation( $config, $def, $summary, $target_dir, $annovar_filter_name );
        }

        addAnnovarMafReport($config, $def, $summary, $target_dir, $annovar_filter_name, $gatk_prefix, $gatk_index, $gatk_index_snv)
      }
    }

    if ($def->{perform_IBS}){
      if ($def->{perform_muTect} || $def->{perform_muTect2indel} || $def->{perform_muTect2}) {
        my $ibs_name = $gatk_prefix . getNextIndex($gatk_index, $gatk_index_snv) . "_IBS";
        $config->{$ibs_name} = {
          class                 => "CQS::ProgramWrapper",
          perform               => 1,
          target_dir            => "${target_dir}/${ibs_name}",
          option                => "",
          interpretor           => "python3",
          program               => "../Variants/IBS.py",
          check_program         => 1,
          docker_prefix         => "gatk4_",
          parameterFile1_arg    => "-i",
          parameterFile1_ref    => $filter_name,
          parameterSampleFile1_arg    => "-f",
          parameterSampleFile1_ref    => "groups",
          output_to_same_folder => 1,
          output_arg            => "-o",
          output_file           => "",
          output_file_ext       => ".ibs_score.mean.csv",
          output_other_ext      => ".ibs_score.csv",
          sh_direct             => 1,
          pbs                   => {
            "email"     => $def->{email},
            "emailType" => $def->{emailType},
            "nodes"     => "1:ppn=1",
            "walltime"  => "10",
            "mem"       => "10gb"
          },
        };
        push @$summary, $ibs_name;
      }
    }

    if ( $def->{perform_vep} ) {
      my $vep_name = $gatk_prefix . getNextIndex($gatk_index, $gatk_index_snv) . "_vep";
      $config->{$vep_name} = {
        class      => "Annotation::Vcf2Maf",
        perform    => 1,
        target_dir => "${target_dir}/$vep_name",
        option     => "",
        source_ref => [ $filter_name, ".vcf" ],
        vcf2maf_pl => getValue( $def, "vcf2maf_pl" ),
        vep_path   => getValue( $def, "vep_path" ),
        vep_data   => getValue( $def, "vep_data" ),
        species    => $species,
        ncbi_build => getValue( $def, "ncbi_build" ),
        filter_vcf => $def->{"vep_filter_vcf"},
        ref_fasta  => $fasta,
        sh_direct  => 1,
        pbs        => {
          "email"    => $email,
          "nodes"    => "1:ppn=" . $max_thread,
          "walltime" => "24",
          "mem"      => "40gb"
        },
      };
      push @$summary, $vep_name;
    }
  }

  if ( $def->{"perform_muTect"} ) {
    my $mutect_index_dic = {};
    my $mutect_index_key = "mutect_Index";
    my $mutect_prefix = "${bam_input}_muTect_";

    my $intervals = $def->{intervals};
    if(!defined $intervals){
      $intervals = $def->{target_intervals_file};
    }

    my $mutectName = $mutect_prefix . getNextIndex($mutect_index_dic, $mutect_index_key) . "_call";
    #print($mutectName);
    $config->{$mutectName} = {
      class        => "GATK::MuTect",
      perform      => 1,
      init_command => $def->{muTect_init_command},
      target_dir   => "${target_dir}/$mutectName",
      option       => getValue( $def, "muTect_option" ),
      source_ref   => [ $bam_input, ".bam\$" ],
      groups_ref   => "groups",
      fasta_file   => $fasta,
      dbsnp_file   => $def->{"dbsnp"},
      intervals    => $intervals,
      bychromosome => 0,
      sh_direct    => 0,
      muTect_jar   => getValue( $def, "muTect_jar" ),
      pbs          => {
        "nodes"    => "1:ppn=1",
        "walltime" => getValue($def, "mutect_walltime", "24"),
        "mem"      => getValue($def, "mutect_memory", "40gb"),
      },
    };
    push @$summary, "${mutectName}";

    add_unique_r($config, $def, $summary, $target_dir, $mutectName . "_qc", "../QC/QCUtils.r,../QC/VcfChromosome.r", ".chromosome.png", [[$mutectName, ".pass.vcf.chromosome.count"]] );

    my $combineVariantsName = $mutect_prefix . getNextIndex($mutect_index_dic, $mutect_index_key) . "_merge";
    $config->{$combineVariantsName} = {
      class                 => "CQS::ProgramWrapper",
      perform               => 1,
      target_dir            => "${target_dir}/${combineVariantsName}",
      option                => "",
      interpretor           => "python3",
      program               => "../GATK/mergeMutect.py",
      check_program         => 1,
      parameterSampleFile1_arg    => "-i",
      parameterSampleFile1_ref    => [ $mutectName, ".pass.vcf.gz\$" ],
      parameterSampleFile1_fileonly  => 0,
      output_to_same_folder => 1,
      output_arg            => "-o",
      output_file_ext       => "_pass.combined.vcf",
      sh_direct             => 1,
      pbs                   => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "10gb"
      },
    };
    push @$summary, $combineVariantsName;

    my $filterVariantsName = $mutect_prefix . getNextIndex($mutect_index_dic, $mutect_index_key) . "_filterDepth";
    $config->{$filterVariantsName} = {
      class                 => "CQS::ProgramWrapper",
      perform               => 1,
      target_dir            => "${target_dir}/${filterVariantsName}",
      option                => "",
      interpretor           => "python3",
      program               => "../GATK/filterMutect.py",
      check_program         => 1,
      parameterFile1_arg    => "-i",
      parameterFile1_ref    => [ $combineVariantsName ],
      output_to_same_folder => 1,
      output_arg            => "-o",
      output_file_ext       => ".filtered.vcf",
      sh_direct             => 1,
      pbs                   => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "10gb"
      },
    };
    push @$summary, $filterVariantsName;

    if ( $def->{perform_annovar} ) {
      my $annovar_name = addAnnovar( $config, $def, $summary, $target_dir, $filterVariantsName, ".vcf\$", $mutect_prefix, $mutect_index_dic, $mutect_index_key );
      if ( $def->{annovar_param} =~ /exac/ || $def->{annovar_param} =~ /1000g/ || $def->{annovar_param} =~ /gnomad/ ) {
        my $annovar_filter_name = addAnnovarFilter( $config, $def, $summary, $target_dir, $annovar_name, $mutect_prefix, $mutect_index_dic, $mutect_index_key);

        if ( defined $def->{annotation_genes} ) {
          $annovar_filter_geneannotation_name = addAnnovarFilterGeneannotation( $config, $def, $summary, $target_dir, $annovar_filter_name );
        }

        addAnnovarMafReport($config, $def, $summary, $target_dir, $annovar_filter_name, $mutect_prefix, $mutect_index_dic, $mutect_index_key)
      }
    }
  }

  if ( $def->{"perform_muTect2"}) {
    my $mutect2_normal_files="mutect2_normal_files";
    my $mutect2_tumor_files="mutect2_tumor_files";

    my $pon = getValue($def, "panel_of_normals", "");
    if($def->{"perform_mutect2_pon"}){
      $pon = add_muTect2_PON($config, $def, $individual, $summary, $target_dir, $bam_input, "PON_muTect2_");
    }

    my $mutect_index_key = "mutect2_key";
    my $mutect_index_dic = { $mutect_index_key => 2 };
    my $mutect_prefix = "muTect2_";

    my $mutect_ref;
    my $mutect2_option = getValue($def, "muTect2_option", "");

    my $perform_mutect2_by_wdl = getValue($def, "perform_mutect2_by_wdl", 1);
    if($perform_mutect2_by_wdl){
      my ($tumor_files, $normal_files) = init_normal_tumor_files($config, $def, $bam_input);

      # print("tumor_files=" . $tumor_files . "\n");
      # print("normal_files=" . $normal_files . "\n");

      my $is_tumor_only = not defined $normal_files;
      #print("is_tumor_only=" . $is_tumor_only . "\n");

      if ($pon eq ""){
        my $mutect2call = addMutect2Wdl($config, $def, $individual, $target_dir, $mutect_prefix, $mutect2_option, 0, $is_tumor_only, $tumor_files, $normal_files, undef, undef);
        $mutect_ref = [ $mutect2call, '.vcf$' ];
      }elsif(-e $pon){
        my $suffix = ($pon =~ /vcf.gz/) ? ".tbi" : ".idx";
        my $mutect2call = addMutect2Wdl($config, $def, $individual, $target_dir, $mutect_prefix, $mutect2_option, 0, $is_tumor_only, $tumor_files, $normal_files, $pon, $pon . $suffix);
        $mutect_ref = [ $mutect2call, '.vcf$' ];
      }else{
        my $mutect2call = addMutect2Wdl($config, $def, $individual, $target_dir, $mutect_prefix, $mutect2_option, 0, $is_tumor_only, $tumor_files, $normal_files, [$pon, 'vcf.gz$'], [$pon, 'vcf.gz.tbi$']);
        $mutect_ref = [ $mutect2call, '.vcf$' ];
      }
    }else{
      my $mutect2call = $mutect_prefix . "01_call";
      addMutect2($config, $def, $individual, $target_dir, $bam_input, $mutect2call, $mutect2_option, 1, $pon);
      $mutect_ref = [ $mutect2call, '\.filtered.vcf.gz$' ];
    }

    my ($annovarMaf,$annovarMafReport) = add_post_mutect($config, $def, $target_dir, $summary, $mutect_prefix, $mutect_index_dic, $mutect_index_key, $mutect_ref);
    #my $mutect2merge = "${bam_input}_muTect2_02_merge";
    #add_combine_mutect($config, $def, $summary, $target_dir, $mutect2merge, [$mutect2call, ".vcf\$"]);
    
    if ($def->{ncbi_build} eq "GRCh38" && $def->{'perform_mutect2_by_wdl'}) {
      my $mutect2callReport = addFilterMafAndReport($config, $def, $summary, $target_dir, $mutect_ref);
      push @$summary, $mutect2callReport;
   
      if (defined($def->{family_info_file}) & defined($def->{patient_info_feature}) ) { #both family_info_file and patient_info_feature defined. Can run Clonevol analysis and ape PhylogeneticTree 
        #make maf of all patients as common_sites and in GATK Intervals format
        my $common_sites=addMafToIntervals( $config, $def, $target_dir, $summary,$mutect_prefix, $mutect_index_dic, $mutect_index_key ,$annovarMaf);
        my $addCollectAllelicCountsCall = addCollectAllelicCounts($config, $def, $summary, $target_dir, $bam_input,$common_sites);

        my $addPrepareClonalAnalysisTask = addAllelicCountsForClonalAnalysis($config, $def, $summary, $target_dir,$addCollectAllelicCountsCall,$mutect2callReport,"cnv_summaryTable");
        my $addSciCloneAndClonevolTask = addSciCloneAndClonevol($config, $def, $summary, $target_dir,$addPrepareClonalAnalysisTask);
        my $addPyCloneVIAndClonevolTask = addPyCloneVIAndClonevol($config, $def, $summary, $target_dir,$addPrepareClonalAnalysisTask);

        my $addApePhylogeneticTreeTask = addApePhylogeneticTree($config, $def, $summary, $target_dir,$mutect2callReport);

        push @$summary, $common_sites,$addCollectAllelicCountsCall,$addPrepareClonalAnalysisTask,$addSciCloneAndClonevolTask,$addPyCloneVIAndClonevolTask,$addApePhylogeneticTreeTask;
      }
    }
  }

  if ( $def->{"perform_cnv_gatk4_somatic"}) {
    my $somaticCNVtask = addSomaticCNV($config, $def, $summary, $target_dir, $bam_input);
  }
  
  if ( $def->{"perform_muTect2indel"} ) {
    my $mutect2indel = "${bam_input}_muTect2indel";
    my $intervals = $def->{muTect2_intervals};
    if(!defined $intervals){
      $intervals = $def->{target_intervals_file};
    }
    $config->{$mutect2indel} = {
      class        => "GATK4::MuTect2indel",
      perform      => 1,
      option       => getValue( $def, "muTect2_option" ),
      init_command => $def->{muTect2_init_command},
      target_dir   => "${target_dir}/$mutect2indel",
      germline_resource => $def->{germline_resource},
      panel_of_normals => $def->{panel_of_normals},
      intervals => $intervals,
      interval_padding => $def->{muTect2_interval_padding},
      source_ref   => [ $bam_input, ".bam\$" ],
      groups_ref   => "groups",
      fasta_file   => $fasta,
      ERC_mode     => $def->{mutect2_ERC_mode},
      sh_direct    => 0,
      pbs          => {
        "nodes"    => "1:ppn=1",
        "walltime" => getValue($def, "mutect2_walltime", "24"),
        "mem"      => getValue($def, "mutect2_memory", "40gb"),
      },
    };
    push @$summary, $mutect2indel;

    add_unique_r($config, $def, $summary, $target_dir, $mutect2indel . "_summary", "../QC/QCUtils.r,../QC/VcfChromosome.r", ".chromosome.png", [[$mutect2indel, ".chromosome.count"]] );

    my $mutect2indel_merge = $mutect2indel . "_merge";
    $config->{$mutect2indel_merge} = {
      class                 => "CQS::ProgramWrapper",
      perform               => 1,
      target_dir            => "${target_dir}/${mutect2indel_merge}",
      option                => "--java-options \"-Xmx40g\" CombineGVCFs \\
      -R $fasta \\
      -V __FILE__ \\
      --OUTPUT __OUTPUT__ \\
      --drop-somatic-filtering-annotations \\
      --input-is-somatic
  ",
      interpretor           => "",
      docker_prefix         => "gatk4_",
      program               => "gatk",
      check_program         => 0,
      source_arg            => "",
      source_type           => "array",
      source_ref            => [ $mutect2indel, ".pass.indel." ],
      source_join_delimiter => " \\\n      -V ",
      output_to_same_folder => 1,
      output_arg            => "",
      output_file_prefix    => ".vcf.gz",
      output_file_ext       => ".vcf.gz",
      sh_direct             => 0,
      pbs                   => {
        "nodes"    => "1:ppn=1",
        "walltime" => "48",
        "mem"      => "40gb"
      },
    };
    push @$summary, $mutect2indel_merge;

    if ( $def->{perform_annovar} ) {
      my $annovar_name = addAnnovar( $config, $def, $summary, $target_dir, $mutect2indel_merge, ".vcf.gz\$" );
    }
  }

  if ( $def->{perform_cnv_cnMOPS} ) {
    my $cnmopsName = "${bam_input}_cnMOPS";
    $config->{$cnmopsName} = {
      class       => "CNV::cnMops",
      perform     => 1,
      target_dir  => "${target_dir}/$cnmopsName",
      option      => "",
      source_ref  => [ $bam_input, ".bam\$" ],
      bedfile     => $covered_bed,
      isbamsorted => 1,
      sh_direct   => 1,
      pbs         => {
        "email"    => $email,
        "nodes"    => "1:ppn=" . $max_thread,
        "walltime" => "24",
        "mem"      => "40gb"
      }
    };
    push @$summary, $cnmopsName;
  }

  my $cnvAnnotationGenesPlot = undef;
  if ( $def->{perform_cnv_gatk4_cohort} ) {
    $cnvAnnotationGenesPlot = addGATK4CNVGermlineCohortAnalysis( $config, $def, $target_dir, [ $bam_input, ".bam\$" ], $bam_input, $individual, $summary, $step3, $step4, $step5, $step6 );
  }

  if ( ( defined $annovar_filter_geneannotation_name ) and ( defined $cnvAnnotationGenesPlot ) ) {
    my $oncoPlotTask = "${bam_input}_SNV_CNV_Oncoplot";
    $config->{$oncoPlotTask} = {
      class                      => "CQS::UniqueR",
      perform                    => 1,
      target_dir                 => $target_dir . '/' . $oncoPlotTask,
      rtemplate                  => "../Visualization/SNV_CNV_OncoPrint.r",
      parameterSampleFile1_ref   => [ $annovar_filter_geneannotation_name, ".oncoprint.tsv\$" ],
      parameterFile1_ref         => [ $cnvAnnotationGenesPlot, ".position.txt.slim" ],
      parameterSampleFile2       => $def->{onco_options},
      parameterSampleFile3       => $def->{onco_sample_groups},
      output_to_result_directory => 1,
      output_file                => "parameterSampleFile1",
      output_file_ext            => ".snv_cnv.txt.png",
      output_other_ext           => ".snv_cnv.txt",
      sh_direct                  => 1,
      'pbs'                      => {
        'nodes'    => '1:ppn=1',
        'mem'      => '40gb',
        'walltime' => '10'
      },
    };
    push @$step6, $oncoPlotTask;
  }

  if ( $def->{perform_cnv_xhmm} ) {
    addXHMM( $config, $def, $target_dir, $bam_ref, $individual, $summary, $summary, $summary, $summary, $summary );
  }
  
  #qc
  if ( getValue( $def, "perform_multiqc" ) ) {
    addMultiQC( $config, $def, $summary, $target_dir, $target_dir );
  }

  #pileup
  $config->{"sequencetask"} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_1 => $individual,
      step_2 => $summary,
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "40gb"
    },
  };

  if ( scalar(@$step3) > 0 ) {
    $config->{"sequencetask"}{"source"}{step_3} = $step3;
  }

  if ( scalar(@$step4) > 0 ) {
    $config->{"sequencetask"}{"source"}{step_4} = $step4;
  }

  if ( scalar(@$step5) > 0 ) {
    $config->{"sequencetask"}{"source"}{step_5} = $step5;
  }

  if ( scalar(@$step6) > 0 ) {
    $config->{"sequencetask"}{"source"}{step_6} = $step6;
  }

  return ($config);
}

sub performExomeSeq {
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

sub performExomeSeqTask {
  my ( $def, $task ) = @_;
  my $config = getConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
