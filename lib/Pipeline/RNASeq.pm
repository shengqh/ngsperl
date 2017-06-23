#!/usr/bin/perl
package Pipeline::RNASeq;

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

our %EXPORT_TAGS = ( 'all' => [qw(performRNASeq performRNASeqTask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  initDefaultValue( $def, "perform_rnaseqc",            0 );
  initDefaultValue( $def, "perform_qc3bam",             0 );
  initDefaultValue( $def, "perform_bamplot",            0 );
  initDefaultValue( $def, "perform_call_variants",      0 );
  initDefaultValue( $def, "aligner",                    "star" );
  initDefaultValue( $def, "use_pearson_in_hca",         1 );
  initDefaultValue( $def, "top25cv_in_hca",             0 );
  initDefaultValue( $def, "use_green_red_color_in_hca", 1 );
  initDefaultValue( $def, "output_bam_to_same_folder",  1 );
  initDefaultValue( $def, "max_thread",                 8 );
  return $def;
}

sub getRNASeqConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

  my $cluster = $def->{cluster};
  my $task    = $def->{task_name};

  my $email = $def->{email};
  my $cqstools = $def->{cqstools} or die "Define cqstools at definition first";
  my $aligner_index;
  my $aligner = $def->{aligner};
  if ( $aligner eq "star" ) {
    $aligner_index = $def->{star_index} or die "Define star_index at definition first";
  }
  else {
    $aligner_index = $def->{hisat2_index} or die "Define hisat2_index at definition first";
  }
  my $transcript_gtf = $def->{transcript_gtf} or die "Define transcript_gtf at definition first";
  my $name_map_file = $def->{name_map_file};

  if ( $def->{perform_rnaseqc} ) {
    defined $def->{rnaseqc_jar} or die "Define rnaseqc_jar first!";
    ( -e $def->{rnaseqc_jar} )  or die "rnaseqc_jar not exists " . $def->{rnaseqc_jar};
    defined $def->{fasta_file}  or die "Define fasta_file for rnaseqc first!";
    ( -e $def->{fasta_file} )   or die "fasta_file not exists " . $def->{fasta_file};
  }

  if ( $def->{perform_qc3bam} ) {
    defined $def->{qc3_perl} or die "Define qc3_perl first!";
    ( -e $def->{qc3_perl} )  or die "qc3_perl not exists " . $def->{qc3_perl};
  }

  if ( $def->{perform_bamplot} ) {
    defined $def->{dataset_name} or die "Define dataset_name for bamplot first!";
    if ( not defined $def->{bamplot_gff} ) {
      defined $def->{gene_names} or die "Define gene_names for bamplot first, seperate by blank space!";
      defined $def->{add_chr}    or die "Define add_chr for bamplot first, check your genome sequence!";
    }
  }

  if ( $def->{perform_call_variants} ) {
    defined $def->{fasta_file}       or die "Define fasta_file for calling variants";
    defined $def->{dbsnp}            or die "Define dbsnp for calling variants";
    defined $def->{gatk_jar}         or die "Define gatk_jar for calling variants";
    defined $def->{picard_jar}       or die "Define picard_jar for calling variants";
    defined $def->{annovar_param}    or die "Define annovar_param for calling variants";
    defined $def->{annovar_db}       or die "Define annovar_db for calling variants";
    defined $def->{annovar_buildver} or die "Define annovar_buildver for calling variants";
  }

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir ) = getPreprocessionConfig($def);

  my $target_dir = $def->{target_dir};

  my $groups_ref = defined $def->{groups} ? "groups" : undef;

  my $configAlignment;
  if ( $aligner eq "star" ) {
    my $starFolder = $target_dir . "/" . getNextFolderIndex($def) . "star";
    $configAlignment = {
      "star" => {
        class                     => "Alignment::STAR",
        perform                   => 1,
        target_dir                => $starFolder,
        option                    => "--twopassMode Basic",
        source_ref                => $source_ref,
        genome_dir                => $aligner_index,
        output_sort_by_coordinate => 1,
        output_to_same_folder     => $def->{output_bam_to_same_folder},
        sh_direct                 => 0,
        pbs                       => {
          "email"    => $email,
          "nodes"    => "1:ppn=" . $def->{max_thread},
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      "star_summary" => {
        class      => "Alignment::STARSummary",
        perform    => 1,
        target_dir => $starFolder,
        option     => "",
        source_ref => [ "star", "_Log.final.out" ],
        sh_direct  => 1,
        pbs        => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      }
    };
    $source_ref = [ "star", "_Aligned.sortedByCoord.out.bam\$" ];
    push @$summary, ("${aligner}_summary");
  }
  else {
    $configAlignment = {
      hisat2 => {
        class                 => "Alignment::Hisat2",
        perform               => 1,
        target_dir            => $target_dir . "/" . getNextFolderIndex($def) . "hisat2",
        option                => "",
        source_ref            => $source_ref,
        genome_dir            => $aligner_index,
        output_to_same_folder => $def->{output_bam_to_same_folder},
        sh_direct             => 1,
        pbs                   => {
          "email"    => $email,
          "nodes"    => "1:ppn=" . $def->{max_thread},
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
    };
    $source_ref = [ "hisat2", ".bam\$" ];
  }

  $configAlignment = merge(
    $configAlignment,
    {
      "featurecount" => {
        class      => "Count::FeatureCounts",
        perform    => 1,
        target_dir => $target_dir . "/" . getNextFolderIndex($def) . "featurecount",
        option     => "-g gene_id -t exon",
        source_ref => $source_ref,
        gff_file   => $transcript_gtf,
        ispairend  => 1,
        sh_direct  => 0,
        pbs        => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      "genetable" => {
        class         => "CQS::CQSDatatable",
        perform       => 1,
        target_dir    => $target_dir . "/" . getNextFolderIndex($def) . "genetable",
        option        => "-k 0 -v 6 -e --fillMissingWithZero",
        source_ref    => "featurecount",
        name_map_file => $name_map_file,
        cqs_tools     => $cqstools,
        sh_direct     => 1,
        pbs           => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "10",
          "mem"      => "10gb"
        },
      },
      "genetable_correlation" => {
        class           => "CQS::UniqueR",
        perform         => 1,
        rCode           => "usePearsonInHCA<-" . $def->{use_pearson_in_hca} . "; useGreenRedColorInHCA<-" . $def->{use_green_red_color_in_hca} . "; top25cvInHCA<-" . $def->{top25cv_in_hca} . "; ",
        target_dir      => $target_dir . "/" . getNextFolderIndex($def) . "genetable_correlation",
        rtemplate       => "countTableVisFunctions.R,countTableGroupCorrelation.R",
        output_file     => "parameterSampleFile1",
        output_file_ext => ".Correlation.png",
        parameterSampleFile1_ref => [ "genetable", ".count\$" ],
        parameterSampleFile2_ref => $groups_ref,
        sh_direct                => 1,
        pbs                      => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      }
    }
  );

  $config = merge( $config, $configAlignment );
  push @$individual, ( "${aligner}", "featurecount" );
  push @$summary,    ( "genetable",  "genetable_correlation" );

  if ( defined $def->{pairs} ) {
    addDEseq2( $config, $def, $summary, "genetable", [ "genetable", ".count\$" ], $def->{target_dir}, $def->{DE_min_median_read} );
  }

  if ( $def->{perform_rnaseqc} ) {
    $config->{rnaseqc} = {
      class          => "QC::RNASeQC",
      perform        => 1,
      target_dir     => $target_dir . "/" . getNextFolderIndex($def) . "rnaseqc",
      init_command   => $def->{rnaseqc_init_command},
      option         => "",
      source_ref     => $source_ref,
      jar            => $def->{rnaseqc_jar},
      fasta_file     => $def->{fasta_file},
      rrna_fasta     => $def->{rrna_fasta},
      transcript_gtf => $transcript_gtf,
      pbs            => {
        "email"    => $email,
        "nodes"    => "1:ppn=" . $def->{max_thread},
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    push( @$summary, "rnaseqc" );
  }

  if ( $def->{perform_qc3bam} ) {
    $config->{qc3} = {
      class          => "QC::QC3bam",
      perform        => 1,
      target_dir     => $target_dir . "/" . getNextFolderIndex($def) . "qc3",
      option         => "",
      transcript_gtf => $transcript_gtf,
      qc3_perl       => $def->{qc3_perl},
      source_ref     => $source_ref,
      pbs            => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    push( @$summary, "qc3" );
  }

  if ( $def->{perform_bamplot} ) {
    if ( not defined $def->{bamplot_gff} ) {
      $config->{gene_pos} = {
        class        => "Annotation::PrepareGenePosition",
        perform      => 1,
        target_dir   => $target_dir . "/" . getNextFolderIndex($def) . "gene_pos",
        option       => "",
        dataset_name => $def->{dataset_name},
        gene_names   => $def->{gene_names},
        add_chr      => $def->{add_chr},
        output_gff   => 1,
        pbs          => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      };
      $config->{bamplot} = {
        class              => "Visualization::Bamplot",
        perform            => 1,
        target_dir         => $target_dir . "/" . getNextFolderIndex($def) . "bamplot",
        option             => "-g " . $def->{dataset_name} . " -y uniform -r --save-temp",
        source_ref         => $source_ref,
        gff_file_ref       => "gene_pos",
        is_rainbow_color   => 0,
        is_single_pdf      => 0,
        is_draw_individual => 0,
        groups             => $def->{"plotgroups"},
        colors             => $def->{"colormaps"},
        sh_direct          => 1,
        pbs                => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      };
      push( @$summary, "gene_pos", "bamplot" );
    }
    else {
      $config->{bamplot} = {
        class              => "Visualization::Bamplot",
        perform            => 1,
        target_dir         => $target_dir . "/" . getNextFolderIndex($def) . "bamplot",
        option             => "-g " . $def->{dataset_name} . " -y uniform -r --save-temp",
        source_ref         => $source_ref,
        gff_file           => $def->{bamplot_gff},
        is_rainbow_color   => 0,
        is_single_pdf      => 0,
        is_draw_individual => 0,
        groups             => $def->{"plotgroups"},
        colors             => $def->{"colormaps"},
        sh_direct          => 1,
        pbs                => {
          "email"    => $email,
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      };
      push( @$summary, "bamplot" );
    }
  }

  if ( $def->{perform_call_variants} ) {
    my $fasta  = getValue( $def, "fasta_file" );
    my $dbsnp  = getValue( $def, "dbsnp" );
    my $gatk   = getValue( $def, "gatk_jar" );
    my $picard = getValue( $def, "picard_jar" );

    $config->{refine} = {
      class              => "GATK::RNASeqRefine",
      perform            => 1,
      target_dir         => $target_dir . "/" . getNextFolderIndex($def) . "refine",
      source_ref         => $source_ref,
      option             => "-Xmx40g",
      fasta_file         => $fasta,
      vcf_files          => [$dbsnp],
      gatk_jar           => $gatk,
      picard_jar         => $picard,
      fixMisencodedQuals => 0,
      replace_read_group => 0,
      reorderChromosome  => 0,
      sorted             => 1,
      sh_direct          => 0,
      pbs                => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };

    $config->{refine_hc} = {
      class         => "GATK::HaplotypeCaller",
      perform       => 1,
      target_dir    => $target_dir . "/" . getNextFolderIndex($def) . "refine_hc",
      option        => "-dontUseSoftClippedBases -stand_call_conf 20.0",
      source_ref    => "refine",
      java_option   => "",
      fasta_file    => $fasta,
      gatk_jar      => $gatk,
      extension     => ".vcf",
      by_chromosome => 0,                                                            #since we have the bed file, we cannot use by_chromosome.
      gvcf          => 0,                                                            #http://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq
      sh_direct     => 0,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };

    $config->{refine_hc_filter} = {
      class       => "GATK::VariantFilter",
      perform     => 1,
      target_dir  => $target_dir . "/" . getNextFolderIndex($def) . "refine_hc_filter",
      option      => "",
      gvcf        => 0,
      vqsr_mode   => 0,
      source_ref  => "refine_hc",
      java_option => "",
      fasta_file  => $fasta,
      dbsnp_vcf   => $dbsnp,
      gatk_jar    => $gatk,
      is_rna      => 1,
      sh_direct   => 1,
      pbs         => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };

    $config->{refine_hc_filter_annovar} = {
      class      => "Annotation::Annovar",
      perform    => 1,
      target_dir => $target_dir . "/" . getNextFolderIndex($def) . "refine_hc_filter_annovar",
      source_ref => "refine_hc_filter",
      option     => $def->{annovar_param},
      annovar_db => $def->{annovar_db},
      buildver   => $def->{annovar_buildver},
      sh_direct  => 1,
      isvcf      => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      },
    };
    push( @$individual, "refine",           "refine" );
    push( @$summary,    "refine_hc_filter", "refine_hc_filter_annovar" );

  }

  $config->{sequencetask} = {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => $individual,
      step2 => $summary,
    },
    sh_direct => 0,
    cluster   => $cluster,
    pbs       => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=" . $def->{max_thread},
      "walltime" => $def->{sequencetask_run_time},
      "mem"      => "40gb"
    },
  };

  return ($config);
}

sub performRNASeq {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getRNASeqConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

sub performRNASeqTask {
  my ( $def, $task ) = @_;

  my $config = getRNASeqConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
