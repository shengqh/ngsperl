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

  initDefaultValue( $def, "perform_gatk_callvariants",   0 );
  initDefaultValue( $def, "gatk_callvariants_vqsr_mode", 1 );
  
  initDefaultValue( $def, "filter_variants_by_allele_frequency",            0 );
  initDefaultValue( $def, "filter_variants_by_allele_frequency_percentage", 0.9 );
  initDefaultValue( $def, "filter_variants_by_allele_frequency_maf",        0.3 );

  initDefaultValue( $def, "perform_muTect",       0 );
  initDefaultValue( $def, "perform_muTect2indel", 0 );
  initDefaultValue( $def, "perform_annovar",      0 );
  initDefaultValue( $def, "perform_cnv",          0 );

  if ( $def->{perform_muTect} || $def->{perform_muTect2indel} ) {
    if ( defined $def->{mills} ) {
      initDefaultValue( $def, "indel_realignment", 1 );
      initDefaultValue( $def, "indel_vcf_files",   $def->{mills} );
    }
    else {
      initDefaultValue( $def, "indel_realignment", 0 );
    }
  }

  initDefaultValue( $def, "perform_multiqc", 1 );

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
  my $cqstools   = getValue( $def, "cqstools" );
  my $max_thread = getValue( $def, "max_thread" );

  my $bam_ref;
  my $fasta;
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

  my $perform_cnv = $def->{perform_cnv_cnMOPs} || $def->{perform_cnv_gatk4_cohort};

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
      remove_duplicate         => 1,
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

    if ( $def->{perform_gatk_callvariants} ) {
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
          "walltime" => "72",
          "mem"      => "40gb"
        },
      };
      push @$individual, ($gvcf_name);

      my $filter_name;
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
          cqstools    => $cqstools,
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

      if ( $def->{filter_variants_by_allele_frequency} ) {
        my $maf_filter_name = $filter_name . "_filterMAF";
        $config->{$maf_filter_name} = {
          class                 => "CQS::ProgramWrapper",
          perform               => 1,
          target_dir            => "${target_dir}/${maf_filter_name}",
          option                => "-p " . $def->{"filter_variants_by_allele_frequency_percentage"} . " -f " . $def->{"filter_variants_by_allele_frequency_maf"},
          interpretor           => "python",
          program               => "../Annotation/filterVcf.py",
          parameterFile1_arg            => "-i",
          parameterFile1_ref            => $filter_name,
          output_to_same_folder => 1,
          output_arg            => "-o",
          output_ext            => ".maf_filtered.vcf",
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

        if ( $def->{annovar_param} =~ /exac/ ) {
          my $annovar_filter_name = addAnnovarFilter( $config, $def, $summary, $target_dir, $annovar_name );

          if ( defined $def->{annotation_genes} ) {
            addAnnovarFilterGeneannotation( $config, $def, $summary, $target_dir, $annovar_filter_name );
          }
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
      my $mutectName = "${refine_name}_muTect";
      $config->{$mutectName} = {
        class        => "GATK::MuTect",
        perform      => 1,
        init_command => $def->{muTect_init_command},
        target_dir   => "${target_dir}/$mutectName",
        option       => getValue( $def, "muTect_option" ),
        java_option  => "-Xmx40g",
        source_ref   => [ $refine_name, ".bam\$" ],
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

      if ( $def->{perform_annovar} ) {
        my $annovar_name = addAnnovar( $config, $def, $summary, $target_dir, $mutectName, ".pass.vcf\$" );
      }
    }

    if ( $def->{"perform_muTect2indel"} ) {
      my $mutect2Name = "${refine_name}_muTect2indel";
      $config->{$mutect2Name} = {
        class        => "GATK::MuTect2Indel",
        perform      => 1,
        init_command => $def->{muTect2_init_command},
        target_dir   => "${target_dir}/$mutect2Name",
        option       => getValue( $def, "muTect2_option" ),
        java_option  => "-Xmx40g",
        source_ref   => [ $refine_name, ".bam\$" ],
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
      my $cnmopsName = "${refine_name}_cnMOPS";
      $config->{$cnmopsName} = {
        class       => "CNV::cnMops",
        perform     => 1,
        target_dir  => "${target_dir}/$cnmopsName",
        option      => "",
        source_ref  => [ $refine_name, ".bam\$" ],
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

    if ( $def->{perform_cnv_gatk4_cohort} ) {
      addGATK4CNVGermlineCohortAnalysis( $config, $def, $target_dir, [ $refine_name, ".bam\$" ], $individual, $summary, $step3, $step4, $step5, $step6 );
    }

    if ( $def->{perform_cnv_xhmm} ) {
      addXHMM( $config, $def, $target_dir, [ $refine_name, ".bam\$" ], $individual, $summary, $step3, $step4, $step5, $step6 );
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
      "walltime" => "72",
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
