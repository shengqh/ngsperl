#!/usr/bin/perl
package Pipeline::MixRNASeq;

use strict;
use warnings;
use List::Util qw(first);
use File::Basename;
use Storable qw(dclone);
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

our %EXPORT_TAGS = ( 'all' => [qw(initializeRNASeqDefaultOptions performMixRNASeq performMixRNASeqTask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeRNASeqDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );

  initDefaultValue( $def, "perform_preprocessing", 1 );
  initDefaultValue( $def, "perform_mapping",       1 );
  initDefaultValue( $def, "perform_counting",      1 );
  initDefaultValue( $def, "perform_count_table",   1 );
  initDefaultValue( $def, "perform_correlation",   1 );
  initDefaultValue( $def, "perform_rnaseqc",       0 );
  initDefaultValue( $def, "perform_qc3bam",        0 );
  initDefaultValue( $def, "perform_bamplot",       0 );
  initDefaultValue( $def, "perform_call_variants", 0 );
  initDefaultValue( $def, "perform_multiqc",       0 );
  initDefaultValue( $def, "perform_webgestalt",    0 );
  initDefaultValue( $def, "perform_webgestaltHeatmap",    0 );
  initDefaultValue( $def, "perform_report",        1 );

  if ( not $def->{"perform_gsea"} ) {
    $def->{"perform_gsea"} = 0;
  }
  else {
    initDefaultValue( $def, "perform_gsea", 0 );
  }

  initDefaultValue( $def, "perform_cutadapt", 0 );

  initDefaultValue( $def, "featureCount_option",        "-g gene_id -t exon" );
  initDefaultValue( $def, "aligner",                    "star" );
  initDefaultValue( $def, "star_option",                "--twopassMode Basic --outSAMmapqUnique 60 --outSAMprimaryFlag AllBestScore" );
  initDefaultValue( $def, "use_pearson_in_hca",         1 );
  initDefaultValue( $def, "top25cv_in_hca",             0 );
  initDefaultValue( $def, "use_green_red_color_in_hca", 1 );
  initDefaultValue( $def, "output_bam_to_same_folder",  1 );
  initDefaultValue( $def, "show_label_PCA",             1 );

  initDefaultValue( $def, "max_thread",            8 );
  initDefaultValue( $def, "sequencetask_run_time", '24' );

  initDefaultValue( $def, "perform_keggprofile",      0 );
  initDefaultValue( $def, "keggprofile_useRawPValue", 0 );
  initDefaultValue( $def, "keggprofile_species",      "hsa" );
  initDefaultValue( $def, "keggprofile_pCut",         0.1 );

  initDefaultValue( $def, "is_paired_end",                   1 );
  initDefaultValue( $def, "DE_pvalue",                       0.05 );
  initDefaultValue( $def, "DE_use_raw_pvalue",               0 );
  initDefaultValue( $def, "DE_fold_change",                  2 );
  initDefaultValue( $def, "DE_export_significant_gene_name", 1 );
  initDefaultValue( $def, "DE_show_gene_cluster",            0 );
  initDefaultValue( $def, "DE_add_count_one",                0 );
  initDefaultValue( $def, "DE_top25only",                    0 );
  initDefaultValue( $def, "DE_detected_in_both_group",       0 );
  initDefaultValue( $def, "DE_perform_wilcox",               0 );
  initDefaultValue( $def, "DE_text_size",                    10 );
  initDefaultValue( $def, "DE_min_median_read",              5 );
  initDefaultValue( $def, "DE_cooksCutoff",                  0.99 );
  initDefaultValue( $def, "perform_DE_proteincoding_gene",   1 );
  initDefaultValue( $def, "perform_proteincoding_gene",      getValue( $def, "perform_DE_proteincoding_gene" ) );

  initDefaultValue( $def, "DE_outputPdf",  getValue( $def, "outputPdf",  0 ) );
  initDefaultValue( $def, "DE_outputPng",  getValue( $def, "outputPng",  1 ) );
  initDefaultValue( $def, "DE_outputTIFF", getValue( $def, "outputTIFF", 0 ) );
  initDefaultValue( $def, "DE_showVolcanoLegend", 1 );

  return $def;
}

sub getDeseq2Suffix {
  my ($config, $def, $deseq2taskname) = @_;

  my $suffix = "";
  if ( ( defined $deseq2taskname ) && ( defined $config->{$deseq2taskname} ) ) {
    if ( getValue( $def, "DE_top25only", 0 ) ) {
      $suffix = $suffix . "_top25";
    }

    if ( getValue( $def, "DE_detected_in_both_group", 0 ) ) {
      $suffix = $suffix . "_detectedInBothGroup";
    }

    my $minMedianInGroup = getValue( $def, "DE_min_median_read", 5 );
    if ( $minMedianInGroup > 0 ) {
      $suffix = $suffix . "_min" . $minMedianInGroup;
    }
    if ( getValue( $def, "DE_use_raw_pvalue", 0 ) ) {
      $suffix = $suffix . "_pvalue" . $def->{DE_pvalue};
    }
    else {
      $suffix = $suffix . "_fdr" . $def->{DE_pvalue};
    }
  }

  return($suffix);
}

sub getMixRNASeqConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeRNASeqDefaultOptions($def);

  my $taskName = $def->{task_name};

  my $email = $def->{email};

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

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  #merge summary and individual 
  push @$individual, @$summary;
  $summary = $individual;

  my $target_dir      = $def->{target_dir};
  my $groups_ref      = defined $def->{groups} ? "groups" : undef;

  $config->{"star"} = {
    class                     => "Alignment::STAR",
    perform                   => 1,
    target_dir                => $target_dir . "/" . getNextFolderIndex($def) . "star",
    option                    => getValue($def, "star_option", ""),
    source_ref                => $source_ref,
    genome_dir                => getValue($def, "star_index"),
    output_sort_by_coordinate => 0,
    output_to_same_folder     => $def->{output_bam_to_same_folder},
    sh_direct                 => 0,
    star_location             => $def->{star_location},
    pbs                       => {
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => "23",
      "mem"       => "40gb"
    },
  };

  $config->{"star_summary"} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . "star_summary",
    option                   => "",
    rtemplate                => "../Alignment/STARFeatureCount.r",
    output_file_ext          => ".STARSummary.csv;.STARSummary.csv.png",
    parameterSampleFile1_ref => [ "star", "_Log.final.out" ],
    sh_direct                => 1,
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    },
  };

  my $filterMixBam_folder = $target_dir . "/" . getNextFolderIndex($def) . "filterMixBam";
  $config->{"filterMixBam"} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => $filterMixBam_folder,
    option                => "--host_prefix " . getValue($def, "host_prefix"),
    interpretor           => "python3",
    check_program         => 1,
    program               => "../Alignment/filterMixBam.py",
    source_ref            => ["star"],
    source_arg            => "-i",
    source_join_delimiter => "",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_file_prefix    => ".fixed.bam",
    output_file_ext       => ".fixed.bam",
    output_other_ext      => ".fixed.txt",
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };

  $config->{"filterMixBam_summary"} = {
    class                     => "CQS::CQSDatatable",
    perform                   => 1,
    target_dir                => $filterMixBam_folder,
    option                    => "-k 0 -v 1",
    source_ref                => ["filterMixBam", ".txt"],
    output_proteincoding_gene => 0,
    sh_direct                 => 1,
    pbs                       => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "40gb"
    },
  };

  $config->{"featureCount"} = {
    class         => "Count::FeatureCounts",
    perform       => 1,
    target_dir    => $target_dir . "/" . getNextFolderIndex($def) . "featureCount",
    option        => "-g gene_id -t exon",
    source_ref    => ["filterMixBam", ".bam\$"],
    gff_file      => getValue($def, "transcript_gtf"),
    is_paired_end => is_paired_end($def),
    sh_direct     => 0,
    pbs           => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "40gb"
    },
  };

  $config->{"featurecount_summary"} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . "featureCountSummary",
    option                   => "",
    rtemplate                => "../Alignment/STARFeatureCount.r",
    output_file_ext          => ".FeatureCountSummary.csv;.FeatureCountSummary.csv.png",
    parameterSampleFile2_ref => [ "featureCount", ".count.summary" ],
    sh_direct                => 1,
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    },
  };

  $config->{"gene_table"} = {
    class                     => "CQS::CQSDatatable",
    perform                   => 1,
    target_dir                => $target_dir . "/" . getNextFolderIndex($def) . "genetable",
    option                    => "-k 0 -v 6 -e --fillMissingWithZero",
    source_ref                => ["featureCount"],
    output_proteincoding_gene => 1,
    name_map_file             => getValue($def, "name_map_file"),
    sh_direct                 => 1,
    pbs                       => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "23",
      "mem"       => "40gb"
    },
  };

  push @$summary, ("star", "star_summary", "filterMixBam", "filterMixBam_summary", "featureCount", "featurecount_summary", "gene_table" );

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      tasks => $summary,
    },
    sh_direct => 0,
    cluster   => $cluster,
    pbs       => {
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => $def->{sequencetask_run_time},
      "mem"       => "40gb"
    },
  };

  return ($config);
}

sub performMixRNASeq {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getMixRNASeqConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

sub performMixRNASeqTask {
  my ( $def, $task ) = @_;

  my $config = getMixRNASeqConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
