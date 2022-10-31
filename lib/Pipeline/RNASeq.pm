#!/usr/bin/perl
package Pipeline::RNASeq;

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

our %EXPORT_TAGS = ( 'all' => [qw(initializeRNASeqDefaultOptions performRNASeq performRNASeqTask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeRNASeqDefaultOptions {
  my $def = shift;

  #print("perform_gsea=" . $def->{perform_gsea} . "\n");

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "ALL" );
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

  #$def->{gsea_jar}        or die "Define gsea_jar at definition first";
  #$def->{gsea_db}         or die "Define gsea_db at definition first";
  #$def->{gsea_categories}

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

sub getRNASeqConfig {
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

  #print(Dumper($def->{groups}));
  #print(Dumper($def->{correlation_groups}));

  my $target_dir      = $def->{target_dir};
  my $groups_ref      = defined $def->{groups} ? "groups" : undef;
  my $aligner         = $def->{aligner};
  my $star_option     = $def->{star_option};
  my $count_table_ref = "files";

  my $multiqc_depedents = $source_ref;

  my $count_file_ref = $def->{count_file};
  if ( $def->{perform_mapping} && $def->{perform_counting} && ( $aligner eq "star" ) && $def->{perform_star_featurecount} ) {
    my $sf_task = addStarFeaturecount($config, $def, $individual, $summary, $target_dir, $source_ref, "" );

    $source_ref      = [ $sf_task, "_Aligned.sortedByCoord.out.bam\$" ];
    $count_table_ref = [ $sf_task, "(?!chromosome).count\$" ];

    $multiqc_depedents = $sf_task;
  }
  else {

    if ( $def->{perform_mapping} ) {
      my $aligner_index;
      if ( $aligner eq "star" ) {
        $aligner_index = $def->{star_index} or die "Define star_index at definition first";
      }
      else {
        $aligner_index = $def->{hisat2_index} or die "Define hisat2_index at definition first";
      }

      my $configAlignment;
      if ( $aligner eq "star" ) {
        my $starFolder = $target_dir . "/" . getNextFolderIndex($def) . "star";
        $configAlignment = {
          "star" => {
            class                     => "Alignment::STAR",
            perform                   => 1,
            target_dir                => $starFolder,
            option                    => $star_option,
            source_ref                => $source_ref,
            genome_dir                => $aligner_index,
            output_sort_by_coordinate => 1,
            output_to_same_folder     => $def->{output_bam_to_same_folder},
            sh_direct                 => 0,
            star_location             => $def->{star_location},
            pbs                       => {
              "nodes"     => "1:ppn=" . $def->{max_thread},
              "walltime"  => "23",
              "mem"       => "40gb"
            },
          },
          "star_summary" => {
            class                    => "CQS::UniqueR",
            perform                  => 1,
            target_dir               => $starFolder,
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
          },
        };

        $source_ref = [ "star", "_Aligned.sortedByCoord.out.bam\$" ];
        push @$summary, ("star_summary");
        $multiqc_depedents = "star";
      }
      else {
        $configAlignment = {
          hisat2 => {
            perform               => 1,
            target_dir            => $target_dir . "/" . getNextFolderIndex($def) . "hisat2",
            class                 => "Alignment::Hisat2",
            option                => "",
            source_ref            => $source_ref,
            genome_dir            => $aligner_index,
            output_to_same_folder => $def->{output_bam_to_same_folder},
            sh_direct             => 1,
            pbs                   => {
              "nodes"     => "1:ppn=" . $def->{max_thread},
              "walltime"  => "23",
              "mem"       => "40gb"
            },
          },
        };
        $source_ref = [ "hisat2", ".bam\$" ];
      }

      $config = merge_hash_right_precedent( $config, $configAlignment );
      push @$individual, $aligner;
      $multiqc_depedents = "hisat2";
    }

    if ( $def->{perform_counting} ) {
      my $transcript_gtf = $def->{transcript_gtf} or die "Define transcript_gtf at definition first";
      if ( $def->{additional_bam_files} ) {
        push @$source_ref, "additional_bam_files";
      }

      my $featureCountFolder = $target_dir . "/" . getNextFolderIndex($def) . "featurecount";
      $config->{"featurecount"} = {
        class         => "Count::FeatureCounts",
        perform       => 1,
        target_dir    => $featureCountFolder,
        option        => "-g gene_id -t exon",
        source_ref    => $source_ref,
        gff_file      => $transcript_gtf,
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
        target_dir               => $featureCountFolder,
        option                   => "",
        rtemplate                => "../Alignment/STARFeatureCount.r",
        #output_file_ext          => ".FeatureCountSummary.csv;.FeatureCountSummary.csv.png",
        output_file_ext          => ".STARSummary.csv;.STARSummary.csv.png",
        parameterSampleFile2_ref => [ "featurecount", ".count.summary" ],
        sh_direct                => 1,
        pbs                      => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "2",
          "mem"       => "10gb"
        },
      };

      push @$individual, "featurecount";
      push @$summary,    "featurecount_summary";
      $count_table_ref = [ "featurecount", ".count\$" ];
      $multiqc_depedents = "featurecount";
    }
  }

  if(getValue($def, "perform_dexseq", 0)){
    my $dexseq_count = "dexseq_count";
    $config->{$dexseq_count} = {
      class        => "Count::DexseqCount",
      perform      => 1,
      target_dir   => $target_dir . "/" . getNextFolderIndex($def) . "$dexseq_count",
      option       => "",
      source_ref   => $source_ref,
      gff_file     => getValue($def, "dexseq_gff"),
      dexseq_count => getValue($def, "dexseq_count.py", "dexseq_count.py"),
      sh_direct    => 0,
      pbs          => {
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    };
    push @$individual, "$dexseq_count";

    my $dexseq_count_table = "dexseq_count_table";
    $config->{$dexseq_count_table} = {
      class         => "CQS::ProgramWrapper",
      perform       => 1,
      target_dir    => $target_dir . "/" . getNextFolderIndex($def) . "$dexseq_count_table",
      interpretor   => "python3",
      program       => "../Count/count_table.py",
      option        => "-p ENS --noheader",
      source_arg    => "-i",
      source_ref    => $dexseq_count,
      output_arg    => "-o",
      output_file_ext => ".count",
      sh_direct     => 1,
      pbs           => {
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "20gb"
      },
    };
    push @$summary, "$dexseq_count_table";
  }

  if(defined $def->{annotation_genes}){
    my $genes_str = $def->{annotation_genes};
    my @genes = split /[;, ]+/, $genes_str;
    my %gene_map = map { $_ => 1 } @genes;
    $config->{annotation_genes} = \%gene_map;
    #print(Dumper($config->{annotation_genes}));

    my $geneLocus = addGeneLocus($config, $def, $summary, $target_dir);

    if($def->{perform_bamsnap}){
      my $bamsnap_task = "annotation_genes_bamsnap";
      addBamsnap($config, $def, $summary, $target_dir, $bamsnap_task, [$geneLocus, "bed"], $source_ref);
    }

    if($def->{bamsnap_coverage}){
      my $coverage_task = "annotation_genes_coverage";
      addGeneCoverage($config, $def, $summary, $target_dir, $coverage_task, "annotation_genes", $source_ref, $geneLocus);
    }

    my $sizeFactorTask = "size_factor";
    addSizeFactor($config, $def, $summary, $target_dir, $sizeFactorTask, $source_ref);
    addPlotGene($config, $def, $summary, $target_dir, "annotation_genes_plot", $sizeFactorTask, [ $geneLocus, ".bed" ], $source_ref);
  }

  if($def->{perform_bamsnap} && $def->{"bamsnap_locus"}){
    addBamsnapLocus($config, $def, $summary, $target_dir, "bamsnap_locus", $source_ref);
  }

  my $perform_count_table = $def->{perform_counting} || $def->{perform_count_table};

  if ($perform_count_table) {
    my $name_map_file = $def->{name_map_file};
    $config->{"genetable"} = {
      class                     => "CQS::CQSDatatable",
      perform                   => 1,
      target_dir                => $target_dir . "/" . getNextFolderIndex($def) . "genetable",
      option                    => "-k 0 -v 6 -e --fillMissingWithZero",
      source_ref                => $count_table_ref,
      output_proteincoding_gene => $def->{perform_proteincoding_gene},
      name_map_file             => $name_map_file,
      sh_direct                 => 1,
      pbs                       => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "23",
        "mem"       => "40gb"
      },
    };

    push @$summary, "genetable";

    if ( getValue( $def, "perform_proteincoding_gene_only", 0 ) ) {
      $count_file_ref = [];
    }
    else {
      $count_file_ref = [ "genetable", "(?<!proteincoding).count\$" ];
    }
    if ( $def->{perform_proteincoding_gene} || getValue( $def, "perform_proteincoding_gene_only", 0 ) ) {
      push @$count_file_ref, "genetable", ".proteincoding.count\$";
    }
  }

  if ( $def->{perform_correlation} ) {
    my $cor_dir = ( defined $config->{genetable} ) ? $config->{genetable}{target_dir} : $target_dir . "/" . getNextFolderIndex($def) . "genetable_correlation";

    my $rCode = getValue( $def, "correlation_rcode", "" );
    $rCode = getOutputFormat( $def, $rCode );
    $rCode = addOutputOption( $def, $rCode, "use_green_red_color_in_hca", $def->{use_green_red_color_in_hca}, "useGreenRedColorInHCA" );
    $rCode = addOutputOption( $def, $rCode, "top25cv_in_hca",             $def->{top25cv_in_hca},             "top25cvInHCA" );

    $config->{"genetable_correlation"} = {
      class           => "CQS::CountTableGroupCorrelation",
      perform         => 1,
      rCode           => $rCode,
      target_dir      => $cor_dir,
      parameterSampleFile4 => {
        "draw_all_groups_in_HCA" => getValue($def, "draw_all_groups_in_HCA", 0),
        "draw_umap" => getValue($def, "draw_umap", 0),
        "heatmap_cexCol" => $def->{heatmap_cexCol},
      },
      rtemplate       => "countTableVisFunctions.R,countTableGroupCorrelation.R",
      output_file     => "parameterSampleFile1",
      output_file_ext => ".Correlation.png;.density.png;.heatmap.png;.PCA.png;.Correlation.Cluster.png",
      sh_direct       => 1,
      pbs             => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "23",
        "mem"       => "40gb"
      },
    };
    if ( is_array($count_file_ref) ) {
      $config->{genetable_correlation}{parameterSampleFile1_ref} = $count_file_ref;
    }
    else {
      $config->{genetable_correlation}{parameterSampleFile1} = { $taskName => [$count_file_ref] };
      $config->{genetable_correlation}{output_to_result_dir} = 1;
    }

    if ( defined $def->{groups} ) {
      $config->{genetable_correlation}{parameterSampleFile2} = { all => $def->{groups} };
    }

    if ( defined $def->{correlation_groups} ) {
      $config->{genetable_correlation}{parameterSampleFile2} = $def->{correlation_groups};
    }

    if ( defined $def->{groups_colors} ) {
      $config->{genetable_correlation}{parameterSampleFile3} = $def->{groups_colors};
    }

    if ( defined $def->{correlation_groups_colors} ) {
      $config->{genetable_correlation}{parameterSampleFile3} = $def->{correlation_groups_colors};
    }

    if ( defined $config->{genetable_correlation}{parameterSampleFile3} ) {
      my $colorGroups = $config->{genetable_correlation}{parameterSampleFile3};
      my $corGroups   = $config->{genetable_correlation}{parameterSampleFile2};
      for my $title ( keys %$corGroups ) {
        my $titleGroups = $corGroups->{$title};
        for my $subGroup ( keys %$titleGroups ) {
          if ( !defined $colorGroups->{$subGroup} ) {
            my %cgroups = %$colorGroups;
            die "Color of group '$subGroup' was not define in " . Dumper($colorGroups);
          }
        }
      }
    }
    push @$summary, "genetable_correlation";

    my $gene_file = $def->{correlation_gene_file};
    if ( defined $gene_file ) {
      $rCode                                                   = $rCode . "suffix<-\"_genes\"; hasRowNames=TRUE;";
      $config->{"genetable_correlation_genes"}                 = dclone( $config->{"genetable_correlation"} );
      $config->{"genetable_correlation_genes"}{target_dir}     = $target_dir . "/" . getNextFolderIndex($def) . "genetable_correlation_genes";
      $config->{"genetable_correlation_genes"}{parameterFile1} = $gene_file;
      $config->{"genetable_correlation_genes"}{rCode}          = $rCode;
      push @$summary, "genetable_correlation_genes";
    }
  }

  my $deseq2taskname;
  my $webgestaltTaskName;
  my $gseaTaskName;
  my $linkTaskName;
  my $webgestaltHeatmapTaskName;
  if ( defined $def->{pairs} ) {
    my $de_prefix;
    my $de_source=$count_file_ref;
    if ( $def->{perform_proteincoding_gene} ) {
      $de_prefix="proteincoding_genetable";
      if(defined $config->{gene_table}){
        $de_source=[ "genetable", ".proteincoding.count\$" ];
      }
    }
    else {
      $de_prefix="genetable";
    }
    $deseq2taskname = addDEseq2( $config, $def, $summary, $de_prefix, $de_source, $def->{target_dir}, $def->{DE_min_median_read} );

    if ( getValue( $def, "perform_webgestalt" ) ) {
      $webgestaltTaskName = addWebgestalt($config, $def, $summary, $target_dir, $deseq2taskname, [ $deseq2taskname, "sig_genename.txt\$" ]);

      #if ( defined $def->{perform_link_webgestalt_deseq2} ) {
      $linkTaskName = $webgestaltTaskName . "_link_deseq2";
      $config->{$linkTaskName} = {
        class                      => "CQS::UniqueR",
        perform                    => 1,
        target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $linkTaskName,
        rtemplate                  => "../Annotation/WebGestaltReportFunctions.r;../Annotation/WebGestaltDeseq2.r",
        rReportTemplate            => "../Annotation/WebGestaltDeseq2.rmd",
        output_to_result_directory => 1,
        output_perSample_file      => "parameterSampleFile1",
        output_perSample_file_ext  => ".html;.html.rds",
        parameterSampleFile1_ref   => [ $webgestaltTaskName, ".txt\$" ],
        parameterSampleFile2_ref   => [ $deseq2taskname, "sig.csv\$" ],
        sh_direct                  => 1,
        rCode                      => "",
        pbs                        => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "23",
          "mem"       => "10gb"
        },
      };
      push( @$summary, $linkTaskName );

      if (getValue( $def, "perform_webgestaltHeatmap" )) {
        $webgestaltHeatmapTaskName = $webgestaltTaskName . "_heatmap_deseq2";
        $config->{$webgestaltHeatmapTaskName} = {
          class                      => "CQS::UniqueR",
          perform                    => 1,
          target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $webgestaltHeatmapTaskName,
          rtemplate                  => "../Annotation/WebGestaltReportFunctions.r;../Annotation/WebGestaltHeatmap.r",
          output_to_result_directory => 1,
          output_perSample_file      => "parameterSampleFile1",
          output_perSample_file_ext  => ".heatmap.png",
          parameterSampleFile1_ref   => [ $webgestaltTaskName, ".txt\$" ],
          parameterSampleFile2_ref   => [ $deseq2taskname, "sig.csv\$" ],
          parameterSampleFile3_ref   => [ $deseq2taskname, "vsd.csv\$" ],
          parameterSampleFile4_ref   => [ $deseq2taskname, ".design\$" ],
          sh_direct                  => 1,
          rCode                      => "",
          pbs                        => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "23",
            "mem"       => "10gb"
          },
        };
        push( @$summary, $webgestaltHeatmapTaskName );
      }
      #}
    }

    #print(Dumper($def));

    if ( getValue( $def, "perform_gsea" ) ) {
      $gseaTaskName = $deseq2taskname . "_GSEA";

      my $pairs = $config->{pairs};
      my $keys = [keys %$pairs];
      #my $suffix = getDeseq2Suffix($config, $def, $deseq2taskname);

      add_gsea($config, $def, $summary, $target_dir, $gseaTaskName, [ $deseq2taskname, "_GSEA.rnk\$" ], $keys, "" );
    }

    if ( $def->{perform_keggprofile} ) {
      my $keggprofile_useRawPValue;
      if ( defined( $def->{keggprofile_useRawPValue} ) ) {
        $keggprofile_useRawPValue = $def->{keggprofile_useRawPValue};
      }
      else {
        die "Define keggprofile_useRawPValue at definition first";
      }
      my $keggprofile_species;
      if ( defined( $def->{keggprofile_species} ) ) {
        $keggprofile_species = $def->{keggprofile_species};
      }
      else {
        die "Define keggprofile_species at definition first";
      }
      my $keggprofile_pCut;
      if ( defined( $def->{keggprofile_pCut} ) ) {
        $keggprofile_pCut = $def->{keggprofile_pCut};
      }
      else {
        die "Define keggprofile_pCut at definition first";
      }
      $config->{keggprofile} = {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $target_dir . "/" . getNextFolderIndex($def) . "keggprofile",
        rtemplate                => "reportFunctions.R;KEGGprofilePerform.R",
        rReportTemplate          => "KEGGprofileReport.Rmd",
        output_file              => "",
        output_file_ext          => ".KEGGprofile.RData",
        parameterSampleFile1_ref => [ $deseq2taskname, "_DESeq2.csv\$" ],
        sh_direct                => 1,
        rCode                    => "useRawPValue=" . $keggprofile_useRawPValue . ";species='" . $keggprofile_species . "';pCut=" . $keggprofile_pCut . ";",
        pbs                      => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "23",
          "mem"       => "10gb"
        },
      };
      push( @$summary, "keggprofile" );
    }
  }

  if ( $def->{perform_rnaseqc} ) {
    my $transcript_gtf = $def->{transcript_gtf} or die "Define transcript_gtf at definition first";
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
        "nodes"     => "1:ppn=" . $def->{max_thread},
        "walltime"  => "23",
        "mem"       => "40gb"
      },
    };
    push( @$summary, "rnaseqc" );
  }

  if ( $def->{perform_rnaseqBamQC} ) {
    my $transcript_gtf = $def->{transcript_gtf} or die "Define transcript_gtf at definition first";

    my $rnaseqBamQC_task = "rnaseqBamQC";
    $config->{$rnaseqBamQC_task} = {
      class                 => "CQS::ProgramWrapperOneToOne",
      perform               => 1,
      target_dir            => "$target_dir/$rnaseqBamQC_task",
      init_command          => "",
      option                => "-i __FILE__ -g '$transcript_gtf' -o __NAME__.txt",
      interpretor           => "python3",
      check_program         => 1,
      program               => "../QC/rnaseqBamQC.py",
      source_ref            => $source_ref,
      source_arg            => "-i",
      output_to_same_folder => 1,
      output_arg            => "-o",
      output_file_prefix    => "",
      output_file_ext       => ".txt",
      sh_direct             => 0,
      pbs                   => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "40gb"
      },
    };
    push (@$individual, $rnaseqBamQC_task);

    my $rnaseqBamQCsummary_task = $rnaseqBamQC_task . "_summary";
    $config->{$rnaseqBamQCsummary_task} = {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $rnaseqBamQCsummary_task,
      rtemplate                => "../QC/rnaseqBamQCsummary.r",
      output_file              => "",
      output_file_ext          => ".txt",
      parameterSampleFile1_ref => $rnaseqBamQC_task,
      sh_direct                => 1,
      rCode                    => "",
      pbs                      => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "2",
        "mem"       => "10gb"
      },
    };
    push( @$summary, $rnaseqBamQCsummary_task );
  }

  if ( $def->{perform_qc3bam} ) {
    my $transcript_gtf = $def->{transcript_gtf} or die "Define transcript_gtf at definition first";
    $config->{qc3} = {
      class          => "QC::QC3bam",
      perform        => 1,
      target_dir     => $target_dir . "/" . getNextFolderIndex($def) . "qc3",
      option         => "",
      transcript_gtf => $transcript_gtf,
      qc3_perl       => $def->{qc3_perl},
      source_ref     => $source_ref,
      pbs            => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "23",
        "mem"       => "40gb"
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
          "nodes"     => "1:ppn=1",
          "walltime"  => "2",
          "mem"       => "10gb"
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
          "nodes"     => "1:ppn=1",
          "walltime"  => "23",
          "mem"       => "10gb"
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
          "nodes"     => "1:ppn=1",
          "walltime"  => "23",
          "mem"       => "10gb"
        },
      };
      push( @$summary, "bamplot" );
    }
  }

  if ( $def->{perform_call_variants} ) {
    my $fasta  = getValue( $def, "fasta_file" );
    my $dbsnp  = getValue( $def, "dbsnp" );

    my $gatk_index = $def;
    my $gatk_index_snv = "SNV_Index";
    my $gatk_prefix;

    my $refine_task;
    my $pass_task;
    if($def->{perform_call_variants_by_gatk4}){
      $refine_task = "gatk4_refine";
      $config->{$refine_task} = {
        class => "CQS::ProgramWrapperOneToOne",
        target_dir => $target_dir . "/" . getNextFolderIndex($def) . $refine_task,
        option => "

gatk MarkDuplicates \\
  --INPUT __FILE__ \\
  --OUTPUT __NAME__.rmdup.bam  \\
  --CREATE_INDEX true \\
  --VALIDATION_STRINGENCY SILENT \\
  --ASSUME_SORTED true \\
  --REMOVE_DUPLICATES true \\
  --METRICS_FILE __NAME__.rmdup.metrics

status=\$?
if [[ \$status -ne 0 ]]; then
  touch __NAME__.rmdup.failed
  rm -f __NAME__.rmdup.bam __NAME__.rmdup.bai
else
  gatk SplitNCigarReads \\
    -R $fasta \\
    -I __NAME__.rmdup.bam \\
    -O __NAME__.rmdup.split.bam 

  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch __NAME__.split.failed
    rm -f __NAME__.rmdup.split.bam __NAME__.rmdup.split.bai
  else
    gatk --java-options \"-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \\
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \\
      -Xloggc:gc_log.log -Xms4000m\" \\
      BaseRecalibrator \\
      -R $fasta \\
      -I __NAME__.rmdup.split.bam \\
      --use-original-qualities \\
      -O __NAME__.rmdup.split.recal.table \\
      -known-sites $dbsnp

    status=\$?
    if [[ \$status -ne 0 ]]; then
      touch __NAME__.recal.failed
      rm -f __NAME__.rmdup.split.recal.table
    else
      gatk --java-options \"-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \\
        -XX:+PrintGCDetails -Xloggc:gc_log.log \\
        -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m\" \\
        ApplyBQSR \\
        --add-output-sam-program-record \\
        -R $fasta \\
        -I __NAME__.rmdup.split.bam \\
        --use-original-qualities \\
        -O __NAME__.rmdup.split.recal.bam \\
        --bqsr-recal-file __NAME__.rmdup.split.recal.table   

      status=\$?
      if [[ \$status -ne 0 ]]; then
        touch __NAME__.recal.failed
        rm -f __NAME__.rmdup.split.recal.bam __NAME__.rmdup.split.recal.bai
      else
        touch __NAME__.recal.succeed
        rm -f __NAME__.rmdup.bam __NAME__.rmdup.bai __NAME__.rmdup.split.bam __NAME__.rmdup.split.bai
        mv __NAME__.rmdup.split.recal.bai __NAME__.rmdup.split.recal.bam.bai
      fi
    fi
  fi
fi

#__OUTPUT__
",
        suffix  => "_rrf",
        docker_prefix => "gatk4_",
        interpretor => "",
        program => "",
        check_program => 0,
        source_arg => "",
        source_ref => $source_ref,
        output_arg => "",
        output_file_ext => ".rmdup.split.recal.bam",
        output_to_same_folder => 1,
        sh_direct   => getValue($def, "sh_direct", 0),
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "23",
          "mem"       => "40gb"
        }
      };

      push( @$individual, $refine_task );

      $gatk_prefix = $refine_task . "_SNV_";

      my $hc_task = $gatk_prefix . getNextIndex($gatk_index, $gatk_index_snv) . "_hc";
      $config->{$hc_task} = {
        class         => "GATK4::HaplotypeCaller",
        perform       => 1,
        target_dir    => $target_dir . "/" . getNextFolderIndex($def) . $hc_task,
        option        => getValue($def, "HaplotypeCaller_option", "--soft-clip-low-quality-ends true --dont-use-soft-clipped-bases true --standard-min-confidence-threshold-for-calling 20"),
        source_ref    => $refine_task,
        java_option   => "",
        bed_file      => $def->{covered_bed},
        fasta_file    => $fasta,
        extension     => ".vcf.gz",
        by_chromosome => 0,                                                            #since we have the bed file, we cannot use by_chromosome.
        gvcf          => 0,                                                            #http://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq
        sh_direct     => 0,
        pbs           => {
          "nodes"     => "1:ppn=8",
          "walltime"  => getValue($def, "HaplotypeCaller_walltime", "23"),
          "mem"       => "40gb"
        },
      };
      push( @$individual, $hc_task );

      my $min_dp = getValue($def, "SNV_minimum_depth", 10);
      my $filter_task = $gatk_prefix . getNextIndex($gatk_index, $gatk_index_snv) . "_filter";
      $config->{$filter_task} = {
        class => "CQS::ProgramWrapperOneToOne",
        target_dir => $target_dir . "/" . $filter_task,
        option => "
gatk VariantFiltration \\
  --R $fasta \\
  --V __FILE__ \\
  --window 35 \\
  --cluster 3 \\
  --filter-name \"FS\" \\
  --filter \"FS > 30.0\" \\
  --filter-name \"DP\" \\
  --filter \"DP < $min_dp\" \\
  --filter-name \"QD\" \\
  --filter \"QD < 2.0\" \\
  -O __NAME__.variant_filtered.vcf.gz

status=\$?
if [[ \$status -ne 0 ]]; then
  touch __NAME__.failed
  rm -f __NAME__.variant_filtered.vcf.gz __NAME__.variant_filtered.vcf.gz.tbi
else
  touch __NAME__.succeed
fi

#__OUTPUT__
",
        suffix  => "_vf",
        docker_prefix => "gatk4_",
        interpretor => "",
        program => "",
        check_program => 0,
        source_arg => "--V",
        source_ref => $hc_task,
        other_localization_ext_array => [".tbi"],
        output_arg => "-O",
        output_file_ext => ".variant_filtered.vcf.gz",
        output_to_same_folder => 1,
        sh_direct   => getValue($def, "sh_direct", 0),
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "10",
          "mem"       => "10gb"
        }
      };
      push( @$individual, $filter_task );

      my $lefttrim_task = $gatk_prefix . getNextIndex($gatk_index, $gatk_index_snv) . "_lefttrim";
      $config->{$lefttrim_task} = {
        class                 => "GATK4::LeftTrim",
        perform               => 1,
        target_dir            => $target_dir . "/" . $lefttrim_task,
        source_ref            => $filter_task,
        option                => "",
        docker_prefix         => "exome_",
        fasta_file            => $fasta,
        extension             => ".variant_filtered.norm.nospan.vcf.gz",
        sh_direct             => 0,
        pbs                   => {
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      };
      push( @$individual, $lefttrim_task );

      $pass_task = $gatk_prefix . getNextIndex($gatk_index, $gatk_index_snv) . "_merge";
      $config->{$pass_task} = {
        class => "CQS::ProgramWrapper",
        target_dir => $target_dir . "/" . $pass_task,
        option => "-i __FILE__ -o __NAME__.pass.vcf.gz

#__OUTPUT__
",
        suffix  => "_mg",
        interpretor => "python3",
        docker_prefix => "exome_",
        program => "../GATK4/combineVCFs.py",
        check_program => 1,
        source_arg => "-i",
        source_ref => $lefttrim_task,
        output_arg => "-o",
        output_file_ext => ".pass.vcf.gz",
        sh_direct   => 1,
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "10",
          "mem"       => "10gb"
        }
      };
      push( @$summary, $pass_task );
    }else{
      my $gatk   = getValue( $def, "gatk_jar" );
      my $picard = getValue( $def, "picard_jar" );

      $refine_task = "refine";
      $gatk_prefix = $refine_task . "_SNV_";

      $config->{$refine_task} = {
        class              => "GATK::RNASeqRefine",
        perform            => 1,
        target_dir         => $target_dir . "/" . getNextFolderIndex($def) . $refine_task,
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
          "nodes"     => "1:ppn=8",
          "walltime"  => "23",
          "mem"       => "40gb"
        },
      };
      push( @$individual, $refine_task );

      my $hc_task = $refine_task . "_hc";
      $config->{$hc_task} = {
        class         => "GATK::HaplotypeCaller",
        perform       => 1,
        target_dir    => $target_dir . "/" . getNextFolderIndex($def) . $hc_task,
        option        => getValue($def, "HaplotypeCaller_option", "--soft-clip-low-quality-ends true --dont-use-soft-clipped-bases true --standard-min-confidence-threshold-for-calling 20"),
        source_ref    => $refine_task,
        java_option   => "",
        fasta_file    => $fasta,
        gatk_jar      => $gatk,
        extension     => ".vcf",
        by_chromosome => 0,                                                            #since we have the bed file, we cannot use by_chromosome.
        gvcf          => 0,                                                            #http://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq
        sh_direct     => 0,
        pbs           => {
          "nodes"     => "1:ppn=8",
          "walltime"  => "23",
          "mem"       => "40gb"
        },
      };
      push( @$individual, $hc_task );

      my $filter_task = $hc_task . "_filter";
      $config->{$filter_task} = {
        class       => "GATK::VariantFilter",
        perform     => 1,
        target_dir  => $target_dir . "/" . getNextFolderIndex($def) . $filter_task,
        option      => "",
        gvcf        => 0,
        vqsr_mode   => 0,
        source_ref  => $hc_task,
        java_option => "",
        fasta_file  => $fasta,
        dbsnp_vcf   => $dbsnp,
        gatk_jar    => $gatk,
        is_rna      => 1,
        sh_direct   => 1,
        pbs         => {
          "nodes"     => "1:ppn=8",
          "walltime"  => "23",
          "mem"       => "40gb"
        },
      };
      push( @$summary, $filter_task );

      $pass_task = $filter_task;
    }

    $multiqc_depedents = $refine_task;

    if ( $def->{filter_variants_by_allele_frequency} ) {
      my $maf_filter_task = $gatk_prefix . getNextIndex($gatk_index, $gatk_index_snv) . "_filterMAF";
      add_maf_filter($config, $def, $summary, $target_dir, $maf_filter_task, $pass_task);
      $pass_task = $maf_filter_task;
    }

    if ( $def->{perform_annovar} ) {
      my $annovar_task = addAnnovar( $config, $def, $summary, $target_dir, $pass_task, undef, $gatk_prefix, $gatk_index, $gatk_index_snv );

      if ( $def->{annovar_param} =~ /exac/ || $def->{annovar_param} =~ /1000g/ || $def->{annovar_param} =~ /gnomad/ ) {
        my $annovar_filter_task = addAnnovarFilter( $config, $def, $summary, $target_dir, $annovar_task, $gatk_prefix, $gatk_index, $gatk_index_snv);

        if ( defined $def->{annotation_genes} ) {
          addAnnovarFilterGeneannotation( $config, $def, $summary, $target_dir, $annovar_filter_task );
        }

        addAnnovarMafReport($config, $def, $summary, $target_dir, $annovar_filter_task, $gatk_prefix, $gatk_index, $gatk_index_snv);
      }
    }
  }

  if ( getValue( $def, "perform_multiqc" ) ) {
    addMultiQC( $config, $def, $summary, $target_dir, $target_dir, $multiqc_depedents );
  }

  if ( getValue( $def, "perform_report" ) ) {
    my @report_files = ();
    my @report_names = ();
    my @copy_files   = ();

    my $version_files = get_version_files($config);

    if ( defined $config->{fastqc_raw_summary} ) {
      push( @report_files, "fastqc_raw_summary",                   ".FastQC.baseQuality.tsv.png" );
      push( @report_files, "fastqc_raw_summary",                   ".FastQC.sequenceGC.tsv.png" );
      push( @report_files, "fastqc_raw_summary",                   ".FastQC.adapter.tsv.png" );
      push( @report_names, "fastqc_raw_per_base_sequence_quality", "fastqc_raw_per_sequence_gc_content", "fastqc_raw_adapter_content" );
    }

    if ( defined $config->{fastqc_post_trim_summary} ) {
      push( @report_files, "fastqc_post_trim_summary",                   ".FastQC.baseQuality.tsv.png" );
      push( @report_files, "fastqc_post_trim_summary",                   ".FastQC.sequenceGC.tsv.png" );
      push( @report_files, "fastqc_post_trim_summary",                   ".FastQC.adapter.tsv.png" );
      push( @report_names, "fastqc_post_trim_per_base_sequence_quality", "fastqc_post_trim_per_sequence_gc_content", "fastqc_post_trim_adapter_content" );
    }

    if ( defined $config->{star_featurecount_summary} ) {
      push( @report_files, "star_featurecount_summary", ".STARSummary.csv.png" );
      push( @report_files, "star_featurecount_summary", ".STARSummary.csv\$" );
      push( @report_names, "STAR_summary",              "STAR_summary_table" );

      push( @report_files, "star_featurecount_summary", ".FeatureCountSummary.csv.png\$" );
      push( @report_files, "star_featurecount_summary", ".FeatureCountSummary.csv\$" );
      push( @report_names, "featureCounts_table_png",   "featureCounts_table" );
    }

    if ( defined $config->{star_summary} ) {
      push( @report_files, "star_summary", ".STARSummary.csv.png" );
      push( @report_files, "star_summary", ".STARSummary.csv\$" );
      push( @report_names, "STAR_summary", "STAR_summary_table" );
    }

    if ( defined $config->{featurecount_summary} ) {
      push( @report_files, "featurecount_summary",    ".FeatureCountSummary.csv.png\$" );
      push( @report_files, "featurecount_summary",    ".FeatureCountSummary.csv\$" );
      push( @report_names, "featureCounts_table_png", "featureCounts_table" );
    }

    if ( defined $config->{genetable} ) {
      push( @copy_files, "genetable", ".count\$", "genetable", ".fpkm.tsv" );
      if($def->{perform_proteincoding_gene}){
        push( @report_files, "genetable",    ".fpkm.proteincoding.tsv\$" );
      }else{
        push( @report_files, "genetable",    ".fpkm.tsv\$" );
      }
      push( @report_names, "genetable_fpkm" );
    }

    if ( defined $config->{genetable_correlation} ) {
      my $suffix = $config->{genetable_correlation}{suffix};
      if ( !defined $suffix ) {
        $suffix = "";
      }
      my $pcoding = $def->{perform_proteincoding_gene} ? ".proteincoding.count" : "";

      my $titles = { "all" => "" };
      if ( is_not_array($count_file_ref) ) {    #count file directly
        $titles->{all} = basename($count_file_ref);
        $pcoding = "";
      }
      if ( defined $config->{genetable_correlation}{parameterSampleFile2} ) {
        my $correlationGroups = $config->{genetable_correlation}{parameterSampleFile2};
        for my $correlationTitle ( keys %$correlationGroups ) {
          my $groups = $correlationGroups->{$correlationTitle};
          if ( is_hash($groups) ) {
            if ( $correlationTitle ne "all" ) {
              $correlationTitle =~ s/\\s+/_/g;
              $titles->{$correlationTitle} = "." . $correlationTitle;
            }
          }
        }
      }

      for my $title ( keys %$titles ) {
        push( @report_files,
          "genetable_correlation", $pcoding . $suffix . $titles->{$title} . ".density.png", "genetable_correlation", $pcoding . $suffix . $titles->{$title} . ".heatmap.png",
          "genetable_correlation", $pcoding . $suffix . $titles->{$title} . ".PCA.png",     "genetable_correlation", $pcoding . $suffix . $titles->{$title} . ".Correlation.Cluster.png" );
        push( @report_names, $title . "_correlation_density", $title . "_correlation_heatmap", $title . "_correlation_PCA", $title . "_correlation_cluster" );
      }
    }

    if ( ( defined $deseq2taskname ) && ( defined $config->{$deseq2taskname} ) ) {
      my $suffix = getDeseq2Suffix($config, $def, $deseq2taskname);      

      my $pairs = $config->{pairs};

      if ( scalar( keys %$pairs ) > 1 ) {
        push( @report_files, $deseq2taskname, "/" . $taskName . ".define.*DESeq2_volcanoPlot.png" );
        push( @report_names, "deseq2_volcano_plot" );
      }
      else {
        push( @report_files, $deseq2taskname, "_DESeq2_volcanoEnhanced.png" );
        push( @report_names, "deseq2_volcano_plot" );
      }
      for my $key ( keys %$pairs ) {
        push( @report_files, $deseq2taskname, "/" . $key . $suffix . "_DESeq2_sig.csv" );
        push( @report_names, "deseq2_" . $key );

        push( @report_files, $deseq2taskname, "/" . $key . ".design" );
        push( @report_names, "deseq2_" . $key . "_design" );

        push( @report_files, $deseq2taskname, "/" . $key . $suffix . "_geneAll_DESeq2-vsd-heatmap.png" );
        push( @report_names, "deseq2_" . $key . "_heatmap" );

        push( @report_files, $deseq2taskname, "/" . $key . $suffix . "_geneAll_DESeq2-vsd-pca.png" );
        push( @report_names, "deseq2_" . $key . "_pca" );
      }
      push( @copy_files, $deseq2taskname, "_DESeq2.csv" );
      push( @copy_files, $deseq2taskname, "_DESeq2_sig.csv" );
      push( @copy_files, $deseq2taskname, "_DESeq2-vsd.csv" );

      #push( @copy_files, $deseq2taskname, "_DESeq2_GSEA.rnk" );
      #push( @copy_files, $deseq2taskname, "_DESeq2_sig_genename.txt" );
      #push( @copy_files, $deseq2taskname, "heatmap.png" );
      #push( @copy_files, $deseq2taskname, "pca.pdf" );
    }

    my $hasFunctionalEnrichment = 0;
    if ( defined $webgestaltTaskName ) {
      push( @copy_files, $webgestaltTaskName, "_geneontology_Biological_Process\$" );
      push( @copy_files, $webgestaltTaskName, "_geneontology_Cellular_Component\$" );
      push( @copy_files, $webgestaltTaskName, "_geneontology_Molecular_Function\$" );
      push( @copy_files, $webgestaltTaskName, "_pathway_KEGG\$" );

      if ( defined $linkTaskName && defined $config->{$linkTaskName} ) {
        push( @copy_files,   $linkTaskName, "txt.html\$" );
      }

      my $pairs = $config->{pairs};
      for my $key ( keys %$pairs ) {
        if ( defined $linkTaskName && defined $config->{$linkTaskName} ) {
          push( @report_files, $linkTaskName, "enrichment_results_" . $key . "_geneontology_Biological_Process.txt.html.rds" );
          push( @report_files, $linkTaskName, "enrichment_results_" . $key . "_geneontology_Cellular_Component.txt.html.rds" );
          push( @report_files, $linkTaskName, "enrichment_results_" . $key . "_geneontology_Molecular_Function.txt.html.rds" );
          push( @report_files, $linkTaskName, "enrichment_results_" . $key . "_pathway_KEGG.txt.html.rds" );
        }
        else {
          push( @report_files, $webgestaltTaskName, "enrichment_results_" . $key . "_geneontology_Biological_Process.txt" );
          push( @report_files, $webgestaltTaskName, "enrichment_results_" . $key . "_geneontology_Cellular_Component.txt" );
          push( @report_files, $webgestaltTaskName, "enrichment_results_" . $key . "_geneontology_Molecular_Function.txt" );
          push( @report_files, $webgestaltTaskName, "enrichment_results_" . $key . "_pathway_KEGG.txt" );
        }
        push( @report_names, "WebGestalt_GO_BP_" . $key );
        push( @report_names, "WebGestalt_GO_CC_" . $key );
        push( @report_names, "WebGestalt_GO_MF_" . $key );
        push( @report_names, "WebGestalt_KEGG_" . $key );
      }
      $hasFunctionalEnrichment = 1;
    }

    if ( defined $gseaTaskName ) {
      push( @copy_files, $gseaTaskName, ".gsea\$" );
      #my $suffix = getDeseq2Suffix($config, $def, $deseq2taskname);      

      my $pairs = $config->{pairs};
      for my $key ( keys %$pairs ) {
        #push( @report_files, $gseaTaskName, "/" . $key . $suffix . "_.*gsea.csv" );
        push( @report_files, $gseaTaskName, "/" . $key . ".gsea.csv" );
        push( @report_names, "gsea_" . $key );
      }
      $hasFunctionalEnrichment = 1;
    }

    my $fcOptions = getValue( $def, "featureCount_option" );
    my $fcMultiMapping = ( $fcOptions =~ /-m/ ) ? "TRUE" : "FALSE";
    my $options = {
      "DE_fold_change"                     => [ getValue( $def, "DE_fold_change",    2 ) ],
      "DE_pvalue"                          => [ getValue( $def, "DE_pvalue",         0.05 ) ],
      "DE_use_raw_pvalue"                  => [ getValue( $def, "DE_use_raw_pvalue", 0 ) ],
      "featureCounts_UseMultiMappingReads" => [$fcMultiMapping],
      "top25cv_in_hca" => [ getValue( $def, "top25cv_in_hca") ? "TRUE" : "FALSE" ],
      "task_name" => $taskName,
      "out.width" => getValue($def, "report.out.width", "60%")
    };

    $config->{report} = {
      class                      => "CQS::BuildReport",
      perform                    => 1,
      target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . "report",
      report_rmd_file            => "../Pipeline/RNASeq.Rmd",
      additional_rmd_files       => "../Pipeline/Pipeline.Rmd;Functions.Rmd",
      parameterSampleFile1_ref   => \@report_files,
      parameterSampleFile1_names => \@report_names,
      parameterSampleFile2       => $options,
      parameterSampleFile3_ref   => \@copy_files,
      parameterSampleFile4       => $version_files,
      parameterSampleFile5       => $def->{software_version},
      parameterSampleFile6       => $def->{groups},
      sh_direct                  => 1,
      pbs                        => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    };
    push( @$summary, "report" );
  }

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
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
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => $def->{sequencetask_run_time},
      "mem"       => "40gb"
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
