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
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performExomeSeq performExomeSeqTask)] );

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

  initDefaultValue( $def, "perform_featureCounts",  0 );

  initDefaultValue( $def, "perform_gatk_callvariants",   0 );
  initDefaultValue( $def, "perform_gatk4_callvariants",  1 );
  initDefaultValue( $def, "gatk_callvariants_vqsr_mode", 1 );
  initDefaultValue( $def, "has_chr_in_chromosome_name" , 0);

  if(defined $def->{"ref_fasta_dict"} && (! defined $def->{"chromosome_names"})){
    my $dictFile = getValue($def, "ref_fasta_dict");
    my $primary_chromosome_only = getValue($def, "primary_chromosome_only", 1);
    my $chromosomes = readChromosomeFromDictFile($dictFile, $primary_chromosome_only);
    initDefaultValue( $def, "chromosome_names" , join(",", @$chromosomes));
  }

  initDefaultValue( $def, "filter_variants_by_allele_frequency",            0 );
  initDefaultValue( $def, "filter_variants_by_allele_frequency_percentage", 0.9 );
  initDefaultValue( $def, "filter_variants_by_allele_frequency_maf",        0.3 );
  initDefaultValue( $def, "filter_variants_fq_equal_1",                     0 );
  
  initDefaultValue( $def, "perform_muTect",       0 );
  initDefaultValue( $def, "perform_muTect2indel", 0 );
  initDefaultValue( $def, "perform_annovar",      0 );
  initDefaultValue( $def, "perform_cnv",          1 );
  initDefaultValue( $def, "perform_vep",          0 );

  if ( $def->{perform_muTect} || $def->{perform_muTect2indel} ) {
    if ( defined $def->{mills} ) {
      initDefaultValue( $def, "indel_realignment", 1 );
      initDefaultValue( $def, "indel_vcf_files",   $def->{mills} );
    }
    else {
      initDefaultValue( $def, "indel_realignment", 0 );
    }
  }

  initDefaultValue( $def, "perform_multiqc", 0 );

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);
  my $step3 = [];
  my $step4 = [];
  my $step5 = [];
  my $step6 = [];

  my $email      = getValue( $def, "email" );
  my $max_thread = getValue( $def, "max_thread" );

  my $geneLocus = undef;
  my $chrCode = getValue($def, "has_chr_in_chromosome_name") ? ";addChr=1" : "";

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
  my $fasta;
  #based on paper https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1097-3, we don't do markduplicate anymore
  if ( $def->{aligner} eq "bwa" ) {
    $fasta = getValue( $def, "bwa_fasta" );
    $config->{ $def->{aligner} } = {
      class                 => "Alignment::BWA",
      perform               => 1,
      target_dir            => "${target_dir}/" . getNextFolderIndex($def) . $def->{aligner},
      option                => getValue( $def, "bwa_option" ),
      bwa_index             => $fasta,
      source_ref            => $source_ref,
      output_to_same_folder => 1,
      picard_jar            => getValue( $def, "picard_jar" ),
      mark_duplicates       => 0,
      sh_direct             => 0,
      pbs                   => {
        "email"    => $email,
        "nodes"    => "1:ppn=" . $max_thread,
        "walltime" => "24",
        "mem"      => "40gb"
      },
    };
    $bam_ref = [ "bwa", ".bam\$" ];
  }
  else {
    die "Unknown alinger " . $def->{aligner};
  }

  push @$individual, ( $def->{aligner} );

  my $perform_cnv = $def->{perform_cnv_cnMOPs} || $def->{perform_cnv_gatk4_cohort} || $def->{perform_cnv_xhmm};

  if ( $def->{perform_gatk_callvariants} || $def->{perform_muTect} || $def->{perform_muTect2_indel} || $perform_cnv ) {
    my $gatk_jar   = getValue( $def, "gatk3_jar" );
    my $picard_jar = getValue( $def, "picard_jar" );

    my $dbsnp = $def->{dbsnp};
    my $indel_vcf_files;
    if ( $def->{indel_realignment} ) {
      if ( !defined $def->{indel_vcf_files} & defined $dbsnp ) {
        $indel_vcf_files = $dbsnp;
      }
      else {
        $indel_vcf_files = getValue( $def, "indel_vcf_files" );
      }
    }
    my $mills = $def->{mills};
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
    my $refine_name = $def->{aligner} . "_refine";
    $config->{$refine_name} = {
      class      => "GATK::Refine",
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
      remove_duplicate         => 0,
      sh_direct                => 0,
      slim_print_reads         => 1,
      samtools_baq_calibration => 0,
      indel_realignment        => $def->{indel_realignment},
      indel_vcf_files          => $indel_vcf_files,
      sorted                   => 1,
      pbs                      => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    };
    push @$individual, ($refine_name);

    my $bam_input = $refine_name;

    if($def->{filter_soft_clip}){
      my $soft_clip_name = $refine_name . "_nosoftclip";
      $config->{$soft_clip_name} = {
        class                 => "CQS::ProgramIndividualWrapper",
        perform               => 1,
        target_dir            => "${target_dir}/${soft_clip_name}",
        option                => "--min-mapq " . getValue($def, "soft_clip_min_mapq", 10),
        interpretor           => "python",
        program               => "../GATK/filterSoftClip.py",
        source_arg            => "-i",
        source_ref            => [ $refine_name, ".bam" ],
        output_to_same_folder => 1,
        output_arg            => "-o",
        output_file_ext       => ".nosoftclip.bam",
        sh_direct             => 0,
        pbs                   => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "10",
          "mem"       => "10gb"
        },
      };
      push @$individual, ($soft_clip_name);
      $bam_input = $soft_clip_name;
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
        output_file_ext          => ".FeatureCountSummary.csv;.FeatureCountSummary.csv.png",
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

    my $filter_name = "";
    if ( $def->{perform_gatk4_callvariants} ) {
      my $gvcf_name         = $bam_input . "_gatk4_hc_gvcf";
      $config->{$gvcf_name} = {
        class             => "GATK4::HaplotypeCaller",
        perform           => 1,
        target_dir        => "${target_dir}/$gvcf_name",
        option            => "",
        source_ref        => $bam_input,
        java_option       => "",
        fasta_file        => $fasta,
        extension         => ".g.vcf",
        bed_file          => $def->{covered_bed},
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

      if(getValue($def, "gatk4_variant_filter_by_chromosome", 0)){
        my $filter_name_chr = $gvcf_name . "_vqsr_chr";
        $config->{$filter_name_chr} = {
          class             => "GATK4::VariantFilterChromosome",
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
            "mem"      => "10gb"
          },
        };
        push @$summary, ($filter_name_chr);

        my $filter_name_chr_recalibrator = $gvcf_name . "_vqsr_chr_recal";
        $config->{$filter_name_chr_recalibrator} = {
          class             => "GATK4::VariantRecalibrator",
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

        my $filter_name_chr_recalibrator_apply = $gvcf_name . "_vqsr_chr_recal_apply";
        $config->{$filter_name_chr_recalibrator_apply} = {
          class             => "GATK4::VariantApplyVQSR",
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

        my $filter_name_chr_recalibrator_apply_gather = $gvcf_name . "_vqsr_chr_recal_apply_gather";
        $config->{$filter_name_chr_recalibrator_apply_gather} = {
          class             => "GATK4::VariantGather",
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
        $filter_name = $gvcf_name . "_vqsr";
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
      my $gvcf_name = $refine_name . "_hc_gvcf";
      $config->{$gvcf_name} = {
        class         => "GATK::HaplotypeCaller",
        perform       => 1,
        target_dir    => "${target_dir}/$gvcf_name",
        option        => "",
        source_ref    => $refine_name,
        java_option   => "",
        fasta_file    => $fasta,
        gatk_jar      => $gatk_jar,
        extension     => ".g.vcf",
        bed_file      => $def->{covered_bed},
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
        $filter_name = $gvcf_name . "_vqsr";
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
        $filter_name = $gvcf_name . "_hardfilter";

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

    my $annovar_filter_geneannotation_name = undef;
    if ( $def->{perform_gatk4_callvariants} or $def->{perform_gatk_callvariants} ) {
      if ( $def->{filter_variants_by_allele_frequency} ) {
        my $maf_filter_name = $filter_name . "_filterMAF";
        $config->{$maf_filter_name} = {
          class                 => "CQS::ProgramWrapper",
          perform               => 1,
          target_dir            => "${target_dir}/${maf_filter_name}",
          option                => "-p " . $def->{"filter_variants_by_allele_frequency_percentage"} . " -f " . $def->{"filter_variants_by_allele_frequency_maf"},
          interpretor           => "python",
          program               => "../Annotation/filterVcf.py",
          parameterFile1_arg    => "-i",
          parameterFile1_ref    => $filter_name,
          output_to_same_folder => 1,
          output_arg            => "-o",
          output_file_ext       => ".maf_filtered.vcf",
          sh_direct             => 1,
          pbs                   => {
            "email"     => $def->{email},
            "emailType" => $def->{emailType},
            "nodes"     => "1:ppn=1",
            "walltime"  => "10",
            "mem"       => "10gb"
          },
        };
        push @$summary, $maf_filter_name;
        $filter_name = $maf_filter_name;
      }

      if ( $def->{perform_annovar} ) {
        my $annovar_name = addAnnovar( $config, $def, $summary, $target_dir, $filter_name, undef );

        if ( $def->{annovar_param} =~ /exac/ || $def->{annovar_param} =~ /1000g/ || $def->{annovar_param} =~ /gnomad/ ) {
          my $annovar_filter_name = addAnnovarFilter( $config, $def, $summary, $target_dir, $annovar_name );

          if ( defined $def->{annotation_genes} ) {
            $annovar_filter_geneannotation_name = addAnnovarFilterGeneannotation( $config, $def, $summary, $target_dir, $annovar_filter_name );
          }

          my $annovar_to_maf = $annovar_filter_name . "_toMAF";
          $config->{$annovar_to_maf} = {
            class      => "Annotation::Annovar2Maf",
            perform    => 1,
            target_dir => $target_dir . "/" . $annovar_to_maf,
            source_ref => [ $annovar_filter_name, "\\.freq0\\..*.filtered.tsv" ],
            refBuild   => getValue( $def, "annovar_buildver" ),
            sh_direct  => 1,
            pbs        => {
              "email"     => $def->{email},
              "emailType" => $def->{emailType},
              "nodes"     => "1:ppn=1",
              "walltime"  => "1",
              "mem"       => "10gb"
            },
          };
          push @$summary, $annovar_to_maf;

          my $annovar_to_maf_report = $annovar_to_maf . "_report";
          $config->{$annovar_to_maf_report} = {
            class                    => "CQS::UniqueR",
            perform                  => 1,
            target_dir               => $target_dir . "/" . $annovar_to_maf_report,
            rtemplate                => "../Annotation/mafReport.r",
            output_file              => "parameterSampleFile1",
            output_file_ext          => ".report.html",
            parameterSampleFile1_ref => [ $annovar_to_maf, ".tsv.maf\$" ],
            parameterFile1           => $def->{family_info_file},
            sh_direct                => 1,
            rCode                    => ( defined $def->{family_info_file} ? "clinicalFeatures=\"" . $def->{family_info_feature} . "\";" : "" ),
#            rCode                    => ( defined $def->{family_info_file} ? "clinicalFeatures=\"" . $def->{family_info_feature} . "\";" : "" )
#              . ( defined $def->{annotation_genes} ? "interestedGeneStr=\"" . $def->{annotation_genes} . "\"" : "" ),
            pbs => {
              "email"     => $def->{email},
              "emailType" => $def->{emailType},
              "nodes"     => "1:ppn=1",
              "walltime"  => "24",
              "mem"       => "10gb"
            },
          };
          push @$summary, $annovar_to_maf_report;
        }
      }

      if ( $def->{perform_vep} ) {
        my $vep_name = $filter_name . "_vep";
        $config->{$vep_name} = {
          class      => "Annotation::Vcf2Maf",
          perform    => 1,
          target_dir => "${target_dir}/$vep_name",
          option     => "",
          source_ref => [ $filter_name, ".vcf" ],
          vcf2maf_pl => getValue( $def, "vcf2maf_pl" ),
          vep_path   => getValue( $def, "vep_path" ),
          vep_data   => getValue( $def, "vep_data" ),
          species    => getValue( $def, "species" ),
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
      my $mutectName = "${bam_input}_muTect";
      $config->{$mutectName} = {
        class        => "GATK::MuTect",
        perform      => 1,
        init_command => $def->{muTect_init_command},
        target_dir   => "${target_dir}/$mutectName",
        option       => getValue( $def, "muTect_option" ),
        java_option  => "-Xmx40g",
        source_ref   => [ $bam_input, ".bam\$" ],
        groups_ref   => "groups",
        fasta_file   => $fasta,
        dbsnp_file   => $def->{"dbsnp"},
        bychromosome => 0,
        sh_direct    => 0,
        muTect_jar   => getValue( $def, "muTect_jar" ),
        pbs          => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "40gb"
        },
      };
      push @$summary, "${mutectName}";

      my $combineVariantsName = $mutectName . "_combined";
      $config->{$combineVariantsName} = {
        class       => "GATK::CombineVariants",
        perform     => 1,
        target_dir  => "${target_dir}/$combineVariantsName",
        option      => "",
        source_ref  => [ $mutectName, ".pass.vcf\$" ],
        java_option => "",
        fasta_file  => $fasta,
        gatk_jar    => $gatk_jar,
        extension   => "_pass.combined.vcf",
        pbs         => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push @$summary, $combineVariantsName;

      if ( $def->{perform_annovar} ) {
        my $annovar_name = addAnnovar( $config, $def, $summary, $target_dir, $combineVariantsName, ".vcf\$" );
        my $annovar_to_maf = $annovar_name . "_toMAF";
        $config->{$annovar_to_maf} = {
          class      => "Annotation::Annovar2Maf",
          perform    => 1,
          target_dir => $target_dir . "/" . $annovar_to_maf,
          source_ref => [$annovar_name],
          refBuild   => getValue( $def, "annovar_buildver" ),
          sh_direct  => 1,
          pbs        => {
            "email"     => $def->{email},
            "emailType" => $def->{emailType},
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        };
        push @$summary, $annovar_to_maf;
      }
    }

    if ( $def->{"perform_muTect2indel"} ) {
      my $mutect2Name = "${bam_input}_muTect2indel";
      $config->{$mutect2Name} = {
        class        => "GATK::MuTect2Indel",
        perform      => 1,
        init_command => $def->{muTect2_init_command},
        target_dir   => "${target_dir}/$mutect2Name",
        option       => getValue( $def, "muTect2_option" ),
        java_option  => "-Xmx40g",
        source_ref   => [ $bam_input, ".bam\$" ],
        groups_ref   => "groups",
        fasta_file   => $fasta,
        dbsnp_file   => $def->{"dbsnp"},
        bychromosome => 0,
        sh_direct    => 0,
        gatk_jar     => getValue( $def, "gatk_jar" ),
        pbs          => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "40gb"
        },
      };
      push @$summary, $mutect2Name;

      if ( $def->{perform_annovar} ) {
        my $annovar_name = addAnnovar( $config, $def, $summary, $target_dir, $mutect2Name, ".pass.vcf\$" );
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
        bedfile     => $def->{covered_bed},
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
        output_to_result_directory => 1,
        output_file                => "parameterSampleFile1",
        output_file_ext            => ".snv_cnv.txt.png;.snv_cnv.txt",
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
      addXHMM( $config, $def, $target_dir, [ $bam_input, ".bam\$" ], $individual, $summary, $step3, $step4, $step5, $step6 );
    }
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
