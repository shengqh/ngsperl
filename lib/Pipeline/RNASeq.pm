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

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );

  initDefaultValue( $def, "perform_preprocessing",          1 );
  initDefaultValue( $def, "perform_mapping",                1 );
  initDefaultValue( $def, "perform_counting",               1 );
  initDefaultValue( $def, "perform_count_table",            1 );
  initDefaultValue( $def, "perform_correlation",            1 );
  initDefaultValue( $def, "perform_rnaseqc",                0 );
  initDefaultValue( $def, "perform_qc3bam",                 0 );
  initDefaultValue( $def, "perform_bamplot",                0 );
  initDefaultValue( $def, "perform_call_variants",          0 );
  initDefaultValue( $def, "perform_call_variants_by_gatk4", 1 );
  initDefaultValue( $def, "perform_multiqc",                0 );
  initDefaultValue( $def, "perform_webgestalt",             0 );
  initDefaultValue( $def, "perform_webgestaltHeatmap",      0 );
  initDefaultValue( $def, "perform_report",                 1 );
  initDefaultValue( $def, "perform_deconvolution",          0 );

  initDefaultValue( $def, "perform_star_fusion", 0 );

  initDefaultValue( $def, "perform_transposable_element", 0 );

  initDefaultValue( $def, "perform_trimmomatic", 0 );
  initDefaultValue( $def, "trimmomatic_option",  ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50" );

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
  if ( $def->{perform_cutadapt} ) {
    #initDefaultValue( $def, "adapter", "CTGTCTCTTATA" );
    initDefaultValue( $def, "min_read_length",                      30 );
    initDefaultValue( $def, "cutadapt_option",                      "-m " . $def->{min_read_length} );
    initDefaultValue( $def, "trim_polyA",                           1 );
    initDefaultValue( $def, "trim_base_quality_after_adapter_trim", 0 );
  } ## end if ( $def->{perform_cutadapt...})

  initDefaultValue( $def, "featureCount_option",        "-g gene_id -t exon" );
  initDefaultValue( $def, "aligner",                    "star" );
  initDefaultValue( $def, "star_option",                "--twopassMode Basic --outSAMmapqUnique 60 --outSAMprimaryFlag AllBestScore" );
  initDefaultValue( $def, "use_pearson_in_hca",         1 );
  initDefaultValue( $def, "top25cv_in_hca",             0 );
  initDefaultValue( $def, "use_green_red_color_in_hca", 0 );
  initDefaultValue( $def, "output_bam_to_same_folder",  1 );

  my $featureCount_option = $def->{featureCount_option};
  if ( $featureCount_option =~ /--fraction/ ) {
    initDefaultValue( $def, "round_count_table", 1 );
  }
  else {
    initDefaultValue( $def, "round_count_table", 0 );
  }

  $def->{"totalCountKey"} = "None";

  if ( defined $def->{files} ) {
    my $files  = $def->{files};
    my $nfiles = keys %$files;
    if ( $nfiles < 20 ) {
      initDefaultValue( $def, "show_label_PCA", 1 );
    }
    else {
      initDefaultValue( $def, "show_label_PCA", 0 );
    }
  } ## end if ( defined $def->{files...})
  else {
    initDefaultValue( $def, "show_label_PCA", 1 );
  }

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
  initDefaultValue( $def, "DE_cooksCutoff",                  "TRUE" );
  initDefaultValue( $def, "DE_independentFiltering",         "TRUE" );
  initDefaultValue( $def, "perform_combatseq",               0 );

  # If we perform_combatseq on count table already, we don't need to perform combatseq in DE analysis.
  if ( $def->{perform_combatseq} ) {
    $def->{"DE_combatseq"} = 0;
  }
  else {
    initDefaultValue( $def, "DE_combatseq", 0 );
  }

  initDefaultValue( $def, "perform_DE_proteincoding_gene", 1 );
  initDefaultValue( $def, "perform_proteincoding_gene",    getValue( $def, "perform_DE_proteincoding_gene" ) );

  initDefaultValue( $def, "DE_outputPdf",         getValue( $def, "outputPdf",  0 ) );
  initDefaultValue( $def, "DE_outputPng",         getValue( $def, "outputPng",  1 ) );
  initDefaultValue( $def, "DE_outputTIFF",        getValue( $def, "outputTIFF", 0 ) );
  initDefaultValue( $def, "DE_showVolcanoLegend", 0 );    #remove size legend from volcano plot

  initDefaultValue( $def, "REPORT_use_enhanced_volcano", 1 );

  initDefaultValue( $def, "minMedian",        1 );
  initDefaultValue( $def, "minMedianInGroup", $def->{DE_min_median_read} );

  return $def;
} ## end sub initializeRNASeqDefaultOptions


sub getDeseq2Suffix {
  my ( $config, $def, $deseq2taskname ) = @_;

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
  } ## end if ( ( defined $deseq2taskname...))

  return ($suffix);
} ## end sub getDeseq2Suffix


sub getRNASeqConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeRNASeqDefaultOptions($def);

  $def->{"perform_link_webgestalt_deseq2_v2"} = 1;

  my $taskName = $def->{task_name};

  my $email = $def->{email};

  if ( $def->{perform_rnaseqc} ) {
    defined $def->{rnaseqc_jar} or die "Define rnaseqc_jar first!";
    ( -e $def->{rnaseqc_jar} )  or die "rnaseqc_jar not exists " . $def->{rnaseqc_jar};
    defined $def->{fasta_file}  or die "Define fasta_file for rnaseqc first!";
    ( -e $def->{fasta_file} )   or die "fasta_file not exists " . $def->{fasta_file};
  } ## end if ( $def->{perform_rnaseqc...})

  if ( $def->{perform_qc3bam} ) {
    defined $def->{qc3_perl} or die "Define qc3_perl first!";
    ( -e $def->{qc3_perl} )  or die "qc3_perl not exists " . $def->{qc3_perl};
  }

  if ( $def->{perform_call_variants} ) {
    defined $def->{fasta_file}       or die "Define fasta_file for calling variants";
    defined $def->{dbsnp}            or die "Define dbsnp for calling variants";
    defined $def->{gatk_jar}         or die "Define gatk_jar for calling variants";
    defined $def->{picard_jar}       or die "Define picard_jar for calling variants";
    defined $def->{annovar_param}    or die "Define annovar_param for calling variants";
    defined $def->{annovar_db}       or die "Define annovar_db for calling variants";
    defined $def->{annovar_buildver} or die "Define annovar_buildver for calling variants";
  } ## end if ( $def->{perform_call_variants...})

  my ( $config, $tasks, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $target_dir = $def->{target_dir};

  if ( $def->{perform_preprocessing_only} ) {

    $config->{sequencetask} = {
      class      => getSequenceTaskClassname($cluster),
      perform    => 1,
      target_dir => "${target_dir}/sequencetask",
      option     => "",
      source     => { tasks => $tasks, },
      sh_direct  => 0,
      cluster    => $cluster,
      pbs        => {
        "nodes"    => "1:ppn=" . $def->{max_thread},
        "walltime" => $def->{sequencetask_run_time},
        "mem"      => "40gb"
      },
    };
    return $config;
  } ## end if ( $def->{perform_preprocessing_only...})

  my $fastq_source_ref = $source_ref;

  #merge summary and individual
  push @$tasks, @$summary;
  $summary = undef;

  #print(Dumper($def->{groups}));
  #print(Dumper($def->{correlation_groups}));

  my $groups_ref      = defined $def->{groups} ? "groups" : undef;
  my $aligner         = $def->{aligner};
  my $star_option     = $def->{star_option};
  my $count_table_ref = "files";

  my $multiqc_depedents = $source_ref;

  if ( $def->{perform_fastq_screen} ) {
    my $fastq_screen_task = "fastq_screen";
    add_fastq_screen( $config, $def, $tasks, $target_dir, $fastq_screen_task, $source_ref );
  }

  my $count_table_column = 6;
  my $count_file_ref     = $def->{count_file};
  if ( $def->{perform_mapping} && $def->{perform_counting} && ( $aligner eq "star" ) && $def->{perform_star_featurecount} ) {
    my $sf_task = addStarFeaturecount( $config, $def, $tasks, $tasks, $target_dir, $source_ref, "" );

    $source_ref      = [ $sf_task, "_Aligned.sortedByCoord.out.bam\$" ];
    $count_table_ref = [ $sf_task, '^(?!.*\.chromosome\.count).*\.count$' ];

    $multiqc_depedents = $sf_task;
  } ## end if ( $def->{perform_mapping...})
  else {
    if ( $def->{perform_mapping} ) {
      if ( $aligner eq "salmon" ) {
        my ( $salmon, $salmon_table ) = add_salmon( $config, $def, $tasks, $tasks, $target_dir, $source_ref, "" );
        $count_table_ref         = [ $salmon, "quant.genes.sf" ];
        $source_ref              = [ $salmon, ".bam\$" ];
        $count_table_column      = 4;
        $def->{perform_counting} = 0;
      } ## end if ( $aligner eq "salmon")
      else {
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
                "nodes"    => "1:ppn=" . $def->{max_thread},
                "walltime" => "23",
                "mem"      => "40gb"
              },
            },
            "star_summary" => {
              class                    => "CQS::UniqueR",
              perform                  => 1,
              target_dir               => $target_dir . "/" . getNextFolderIndex($def) . "star_summary",
              option                   => "",
              rtemplate                => "../Alignment/AlignmentUtils.r,../Alignment/STARFeatureCount.r",
              output_file_ext          => ".STARSummary.csv;.STARSummary.csv.png",
              parameterSampleFile1_ref => [ "star", "_Log.final.out" ],
              sh_direct                => 1,
              pbs                      => {
                "nodes"    => "1:ppn=1",
                "walltime" => "2",
                "mem"      => "10gb"
              },
            },
          };

          $source_ref = [ "star", "_Aligned.sortedByCoord.out.bam\$" ];
          push @$tasks, ("star_summary");
        } ## end if ( $aligner eq "star")
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
                "nodes"    => "1:ppn=" . $def->{max_thread},
                "walltime" => "23",
                "mem"      => "40gb"
              },
            },
          };
          $source_ref = [ "hisat2", ".bam\$" ];
        } ## end else [ if ( $aligner eq "star")]

        $config = merge_hash_right_precedent( $config, $configAlignment );
        push @$tasks, $aligner;

        $multiqc_depedents = $source_ref;
        if ( $def->{perform_umitools} ) {
          my $dedup_task   = "umitools_dedup";
          my $dedup_option = getValue( $def, "unitools_dedup_option", "--method=unique" );

          my $pairend_option = is_paired_end($def) ? "--paired" : "";
          $config->{$dedup_task} = {
            class      => "CQS::ProgramWrapperOneToOne",
            target_dir => $target_dir . "/" . getNextFolderIndex($def) . $dedup_task,
            option     => "
umi_tools dedup $dedup_option $pairend_option --output-stats __NAME__ --stdin __FILE__ --stdout __NAME__.dedup.bam

samtools index __NAME__.dedup.bam

samtools flagstat __NAME__.dedup.bam > __NAME__.dedup.bam.flagstat

#__OUTPUT__
",
            interpretor           => "",
            check_program         => 0,
            program               => "",
            source_arg            => "--stdin",
            source_ref            => $source_ref,
            output_arg            => "--stdout",
            output_file_prefix    => ".dedup.bam",
            output_file_ext       => ".dedup.bam",
            output_to_same_folder => 1,
            docker_prefix         => "umitools_",
            #no_docker => 1,
            use_tmp_folder => 0,
            sh_direct      => 0,
            pbs            => {
              "nodes"    => "1:ppn=1",
              "walltime" => "10",
              "mem"      => "10gb"
            }
          };
          $source_ref = $dedup_task;
          push( @$tasks, $dedup_task );
        } ## end if ( $def->{perform_umitools...})
      } ## end else [ if ( $aligner eq "salmon")]

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
          option        => getValue( $def, "featureCount_option", "-g gene_id -t exon" ),
          source_ref    => $source_ref,
          gff_file      => $transcript_gtf,
          is_paired_end => is_paired_end($def),
          sh_direct     => 0,
          pbs           => {
            "nodes"    => "1:ppn=1",
            "walltime" => "23",
            "mem"      => "40gb"
          },
        };
        $config->{"featurecount_summary"} = {
          class                    => "CQS::UniqueR",
          perform                  => 1,
          target_dir               => "${featureCountFolder}_summary",
          option                   => "",
          rtemplate                => "../Alignment/AlignmentUtils.r,../Alignment/STARFeatureCount.r",
          output_file_ext          => ".FeatureCountSummary.csv;.FeatureCountSummary.csv.png",
          parameterSampleFile2_ref => [ "featurecount", ".count.summary" ],
          sh_direct                => 1,
          pbs                      => {
            "nodes"    => "1:ppn=1",
            "walltime" => "2",
            "mem"      => "10gb"
          },
        };

        push @$tasks, "featurecount";
        push @$tasks, "featurecount_summary";
        $count_table_ref   = [ "featurecount", ".count\$" ];
        $multiqc_depedents = "featurecount";
      } ## end if ( $def->{perform_counting...})
    } ## end if ( $def->{perform_mapping...})
  } ## end else [ if ( $def->{perform_mapping...})]

  if ( $def->{perform_star_fusion} ) {
    my $star_fusion_task = "star_fusion";
    my $genome_lib_dir   = getValue( $def, "star_fusion_genome_lib_dir" );
    my $thread           = $def->{max_thread};
    my $memory_gb        = getValue( $def, "star_fusion_memory_gb", "40" );
    my $memory_sort      = $memory_gb - 5;
    $config->{$star_fusion_task} = {
      class         => "CQS::ProgramWrapperOneToOne",
      perform       => 1,
      target_dir    => $target_dir . "/" . getNextFolderIndex($def) . $star_fusion_task,
      program       => "",
      check_program => 0,
      option        => "
     
STAR-Fusion --genome_lib_dir $genome_lib_dir \\
  --left_fq __FILE__ \\
  --output_dir . \\
  --CPU $thread \\
  --STAR_limitBAMsortRAM ${memory_sort}G \\
  --FusionInspector validate \\
  --examine_coding_effect \\
  --STAR_SortedByCoordinate \\
  --denovo_reconstruct

status=\$?
if [[ \$status -eq 0 ]]; then
  mv star-fusion.fusion_predictions.abridged.coding_effect.tsv __NAME___star-fusion.fusion_predictions.abridged.coding_effect.tsv
  mv star-fusion.fusion_predictions.abridged.tsv __NAME___star-fusion.fusion_predictions.abridged.tsv 
  mv star-fusion.fusion_predictions.tsv __NAME___star-fusion.fusion_predictions.tsv
  touch __NAME__.star_fusion.succeed
  rm -f __NAME__.star_fusion.failed

  if [[ -s FusionInspector-validate/finspector.fusion_inspector_web.html ]]; then
    mv FusionInspector-validate/finspector.fusion_inspector_web.html FusionInspector-validate/__NAME___finspector.fusion_inspector_web.html
    mv FusionInspector-validate/finspector.FusionInspector.fusions.tsv FusionInspector-validate/__NAME___finspector.FusionInspector.fusions.tsv
  fi
else
  touch __NAME__.star_fusion.failed
  rm -f __NAME__.star_fusion.succeed
fi

rm -rf Aligned.sortedByCoord.out.bam Aligned.out.bam Chimeric.out.junction _* star-fusion.preliminary

STAR-Fusion --version | grep version | cut -d ':' -f2 | awk '{print \"STAR-Fusion,v\"\$1}' > __NAME__.version

",
      parameterSampleFile1_join_delimiter => " \\\n  --right_fq ",
      parameterSampleFile1_ref            => $fastq_source_ref,
      output_ext                          => "_star-fusion.fusion_predictions.tsv",
      output_file_ext                     => ".version,_star-fusion.fusion_predictions.abridged.coding_effect.tsv",
      docker_prefix                       => "star_fusion_",
      output_to_same_folder               => 0,
      no_output                           => 1,
      sh_direct                           => 0,
      pbs                                 => {
        "nodes"    => "1:ppn=" . $def->{max_thread},
        "walltime" => "23",
        "mem"      => $memory_gb . "gb"
      },
    };
    push @$tasks, $star_fusion_task;

    my $star_fusion_summary_task = "star_fusion_summary";
    $config->{$star_fusion_summary_task} = {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $star_fusion_summary_task,
      option                   => "",
      rtemplate                => "../Fusion/star_fusion_summary.r",
      output_file_ext          => ".fusion_count.csv,.fusion_count.png",
      parameterSampleFile1_ref => [ $star_fusion_task, ".fusion_predictions.tsv" ],
      sh_direct                => 1,
      pbs                      => {
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      },
    };
    push @$tasks, $star_fusion_summary_task;
  } ## end if ( $def->{perform_star_fusion...})

  if ( getValue( $def, "perform_dexseq", 0 ) ) {
    my $dexseq_count = "dexseq_count";
    $config->{$dexseq_count} = {
      class        => "Count::DexseqCount",
      perform      => 1,
      target_dir   => $target_dir . "/" . getNextFolderIndex($def) . "$dexseq_count",
      option       => "",
      source_ref   => $source_ref,
      gff_file     => getValue( $def, "dexseq_gff" ),
      dexseq_count => getValue( $def, "dexseq_count.py", "dexseq_count.py" ),
      sh_direct    => 0,
      pbs          => {
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    };
    push @$tasks, "$dexseq_count";

    my $dexseq_count_table = "dexseq_count_table";
    $config->{$dexseq_count_table} = {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      target_dir               => $target_dir . "/" . getNextFolderIndex($def) . "$dexseq_count_table",
      rtemplate                => "../Count/dexseq_count_table.r",
      parameterSampleFile1_ref => $dexseq_count,
      parameterFile1           => getValue( $def, "name_map_file" ),
      output_file_ext          => ".exon.count",
      sh_direct                => 1,
      pbs                      => {
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "20gb"
      },
    };
    push @$tasks, "$dexseq_count_table";
  } ## end if ( getValue( $def, "perform_dexseq"...))

  if ( defined $def->{annotation_genes_bed} ) {
    my $genes_file = $def->{annotation_genes_bed};

    my $gene_map = {};
    open( FH, '<', $genes_file ) or die $!;
    while (<FH>) {
      chomp;
      my @parts = split(/\t/);
      $gene_map->{ $parts[4] } = 1;
    }
    close(FH);

    $config->{annotation_genes}     = $gene_map;
    $config->{annotation_genes_ref} = { $taskName => [$genes_file] };

    my $geneLocus = "annotation_genes_ref";

    if ( $def->{perform_bamsnap} ) {
      my $bamsnap_task = "annotation_genes_bamsnap";
      addBamsnap( $config, $def, $tasks, $target_dir, $bamsnap_task, $geneLocus, $source_ref );
    }

    if ( $def->{bamsnap_coverage} ) {
      my $coverage_task = "annotation_genes_coverage";
      addGeneCoverage( $config, $def, $tasks, $target_dir, $coverage_task, "annotation_genes", $source_ref, $geneLocus );
    }

    my $sizeFactorTask = "size_factor";
    addSizeFactor( $config, $def, $tasks, $target_dir, $sizeFactorTask, $source_ref );
    addPlotGene( $config, $def, $tasks, $target_dir, "annotation_genes_plot", $sizeFactorTask, [ $geneLocus, ".bed" ], $source_ref );
  } ## end if ( defined $def->{annotation_genes_bed...})

  if ( defined $def->{annotation_genes} ) {
    my $genes_str = $def->{annotation_genes};
    $genes_str =~ s/\s+/,/g;
    my @genes    = split /[;, ]+/, $genes_str;
    my %gene_map = map { $_ => 1 } @genes;
    $config->{annotation_genes} = \%gene_map;
    #print(Dumper($config->{annotation_genes}));

    my $geneLocus = addGeneLocus( $config, $def, $tasks, $target_dir, "annotation_genes" );

    if ( $def->{perform_bamsnap} ) {
      my $bamsnap_task = "annotation_genes_bamsnap";
      addBamsnap( $config, $def, $tasks, $target_dir, $bamsnap_task, [ $geneLocus, "bed" ], $source_ref );
    }

    if ( $def->{bamsnap_coverage} ) {
      my $coverage_task = "annotation_genes_coverage";
      addGeneCoverage( $config, $def, $tasks, $target_dir, $coverage_task, "annotation_genes", $source_ref, $geneLocus );
    }

    my $sizeFactorTask = "size_factor";
    addSizeFactor( $config, $def, $tasks, $target_dir, $sizeFactorTask, $source_ref );
    addPlotGene( $config, $def, $tasks, $target_dir, "annotation_genes_plot", $sizeFactorTask, [ $geneLocus, ".bed" ], $source_ref );
  } ## end if ( defined $def->{annotation_genes...})

  if ( $def->{perform_bamsnap} && $def->{"bamsnap_locus"} ) {
    addBamsnapLocus( $config, $def, $tasks, $target_dir, "bamsnap_locus", $source_ref );
  }

  my $perform_count_table = $def->{perform_counting} || $def->{perform_count_table};

  if ($perform_count_table) {
    my $name_map_file = $def->{name_map_file};
    $config->{"genetable"} = {
      class                     => "CQS::CQSDatatable",
      perform                   => 1,
      target_dir                => $target_dir . "/" . getNextFolderIndex($def) . "genetable",
      option                    => "-k 0 -v $count_table_column -e --fillMissingWithZero",
      source_ref                => $count_table_ref,
      remove_chrM_genes         => getValue( $def, "remove_chrM_genes", 0 ),
      round_count_table         => getValue( $def, "round_count_table", 0 ),
      output_proteincoding_gene => $def->{perform_proteincoding_gene},
      name_map_file             => $name_map_file,
      sh_direct                 => 1,
      pbs                       => {
        "nodes"    => "1:ppn=1",
        "walltime" => "23",
        "mem"      => "40gb"
      },
    };

    push @$tasks, "genetable";

    if ( getValue( $def, "perform_proteincoding_gene_only", 0 ) ) {
      $count_file_ref = [];
    }
    else {
      $count_file_ref = [ "genetable", "(?<!proteincoding).count\$" ];
    }
    if ( $def->{perform_proteincoding_gene} || getValue( $def, "perform_proteincoding_gene_only", 0 ) ) {
      push @$count_file_ref, "genetable", ".proteincoding.count\$";
    }
  } ## end if ($perform_count_table)

  # perform correlation before possible combatseq
  if ( $def->{perform_correlation} ) {
    my $cor_dir = ( defined $config->{genetable} ) ? $config->{genetable}{target_dir} : $target_dir . "/" . getNextFolderIndex($def) . "genetable_correlation";
    add_table_correlation( $config, $def, $tasks, "genetable_correlation", $cor_dir, $count_file_ref );
  }

  if ( getValue( $def, "perform_combatseq", 0 ) ) {
    my $combatseq_task  = "genetable_combatseq";
    my $output_file_ext = ".combatseq.count,.sva.version";
    if ( scalar(@$count_file_ref) == 4 ) {    #has both count and proteincoding.count
      $output_file_ext = $output_file_ext . ",.combatseq.proteincoding.count";
    }
    $config->{$combatseq_task} = {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $combatseq_task,
      option                   => "",
      rtemplate                => "../Count/combatseq.r",
      output_file_ext          => $output_file_ext,
      parameterSampleFile1_ref => $count_file_ref,
      parameterFile2           => $def->{covariance_file},
      docker_prefix            => "correlation_",
      sh_direct                => 1,
      pbs                      => {
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      },
    };
    push @$tasks, $combatseq_task;

    if ( $def->{perform_correlation} ) {
      my $cor_dir  = $config->{$combatseq_task}{target_dir};
      my $cor_task = $combatseq_task . "_correlation";
      add_table_correlation( $config, $def, $tasks, $cor_task, $cor_dir, [ $combatseq_task, ".count" ] );
    }

    $count_file_ref = [ $combatseq_task, ".count\$" ];
  } ## end if ( getValue( $def, "perform_combatseq"...))

  if ( getValue( $def, "perform_transposable_element", 0 ) ) {
    my $star_index = $def->{star_index} or die "Define star_index at definition first";

    my $transcript_gtf = $def->{transcript_gtf} or die "Define transcript_gtf at definition first";
    my $name_map_file  = $def->{name_map_file}  or die "Define name_map_file at definition first";

    my $te_rmsk_gtf = $def->{transposable_element_rmsk_gtf} or die "Define transposable_element_rmsk_gtf at definition first";

    $config->{"te_star"} = {
      class                     => "Alignment::STAR",
      perform                   => 1,
      target_dir                => $target_dir . "/" . getNextFolderIndex($def) . "te_star",
      option                    => $star_option . " --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100 ",
      source_ref                => $fastq_source_ref,
      genome_dir                => $star_index,
      output_sort_by_coordinate => 1,
      output_to_same_folder     => $def->{output_bam_to_same_folder},
      sh_direct                 => 0,
      star_location             => $def->{star_location},
      pbs                       => {
        "nodes"    => "1:ppn=" . $def->{max_thread},
        "walltime" => "23",
        "mem"      => "40gb"
      },
    };
    $config->{"te_star_summary"} = {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      target_dir               => $target_dir . "/" . getNextFolderIndex($def) . "te_star_summary",
      option                   => "",
      rtemplate                => "../Alignment/AlignmentUtils.r,../Alignment/STARFeatureCount.r",
      output_file_ext          => ".STARSummary.csv;.STARSummary.csv.png",
      parameterSampleFile1_ref => [ "te_star", "_Log.final.out" ],
      sh_direct                => 1,
      pbs                      => {
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      },
    };

    my $te_source_ref = [ "te_star", "_Aligned.sortedByCoord.out.bam\$" ];
    push @$tasks, ( "te_star", "te_star_summary" );

    $config->{"te_star_count"} = {
      class         => "CQS::ProgramWrapperOneToOne",
      perform       => 1,
      target_dir    => $target_dir . "/" . getNextFolderIndex($def) . "te_star_count",
      interpretor   => "",
      program       => "",
      check_program => 0,
      option        => "
TEcount --version | cut -d ' ' -f 2 > __NAME__.version

TEcount \\
  --BAM __FILE__ \\
  --GTF $transcript_gtf \\
  --TE $te_rmsk_gtf \\
  --sortByPos \\
  --project __NAME__
",
      source_arg      => "-i",
      source_ref      => $te_source_ref,
      no_output       => 1,
      output_file_ext => ".cntTable",
      sh_direct       => 0,
      docker_prefix   => "tetranscripts_",
      pbs             => {
        "nodes"    => "1:ppn=8",
        "walltime" => "10",
        "mem"      => "20gb"
      },
    };
    push @$tasks, "te_star_count";

    $config->{"te_genetable"} = {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      target_dir               => $target_dir . "/" . getNextFolderIndex($def) . "te_genetable",
      option                   => "",
      rtemplate                => "../Count/TECountTable.r",
      parameterSampleFile1_ref => [ "te_star_count", ".cntTable" ],
      parameterFile1           => $name_map_file,
      parameterFile2_ref       => $count_file_ref,
      output_file_ext          => ".count",
      sh_direct                => 1,
      pbs                      => {
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      },
    };

    push @$tasks, "te_genetable";

    if ( defined $def->{pairs} ) {
      my $de_prefix      = "te_genetable";
      my $deseq2taskname = addDEseq2(
        $config,
        $def,
        $tasks,
        $de_prefix,
        [ "te_genetable", ".count" ],
        $def->{target_dir},
        $def->{DE_min_median_read},
        undef,    #$libraryFile,
        undef,    #$libraryKey,
        undef,    #$feature_name_regex,
        undef,    #$n_first
        "transposable elements"
      );
    } ## end if ( defined $def->{pairs...})

    if ( getValue( $def, "perform_transposable_element_bamplot", 0 ) ) {
      my $te_rmsk_bed = $def->{transposable_element_rmsk_bed} or die "Define transposable_element_rmsk_bed at definition first";
      my $te_names    = getValue( $def, "transposable_element_bamplot_names" );
      $config->{"te_bamplot_gff"} = {
        class           => "CQS::ProgramWrapperOneToOne",
        perform         => 1,
        target_dir      => $target_dir . "/" . getNextFolderIndex($def) . "te_bamplot_gff",
        interpretor     => "python3",
        program         => "../Format/bed2gff.py",
        check_program   => 1,
        option          => "-b $te_rmsk_bed -n \'$te_names\' -o __NAME__.gff",
        source_arg      => "-b",
        source          => { getValue( $def, "task_name" ) => [$te_rmsk_bed] },
        no_output       => 1,
        output_file_ext => ".gff",
        sh_direct       => 0,
        pbs             => {
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      };
      push @$tasks, "te_bamplot_gff";

      my $plotgroups = $def->{plotgroups};
      if ( !defined $plotgroups ) {
        my $files         = getValue( $config, "files" );
        my @sortedSamples = sort keys %$files;
        $plotgroups = { getValue( $def, "task_name" ) => \@sortedSamples };
      }
      $config->{plotgroups} = $plotgroups;

      $config->{"te_bamplot"} = {
        class              => "Visualization::Bamplot",
        perform            => 1,
        target_dir         => "${target_dir}/" . getNextFolderIndex($def) . "te_bamplot",
        option             => getValue( $def, "bamplot_option" ),
        source_ref         => $te_source_ref,
        groups_ref         => "plotgroups",
        gff_file_ref       => "te_bamplot_gff",
        is_rainbow_color   => 0,
        is_draw_individual => 0,
        is_single_pdf      => 1,
        draw_by_r          => getValue( $def, "bamplot_draw_by_r",        1 ),
        draw_by_r_width    => getValue( $def, "bamplot_draw_by_r_width",  10 ),
        draw_by_r_height   => getValue( $def, "bamplot_draw_by_r_height", 10 ),
        sh_direct          => 1,
        pbs                => {
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      };
      push @$tasks, ("te_bamplot");

    } ## end if ( getValue( $def, "perform_transposable_element_bamplot"...))
  } ## end if ( getValue( $def, "perform_transposable_element"...))

  if ( $def->{perform_deconvolution} ) {
    my $cpm_task         = "deconvolution_1_cpm";
    my $basicMatrix_file = getValue( $def, "deconvolution_basicMatrix" );
    $config->{$cpm_task} = {
      class              => "CQS::UniqueR",
      perform            => 1,
      target_dir         => "$target_dir/$cpm_task",
      option             => "",
      rtemplate          => "../Deconvolution/cpm.r",
      parameterFile1_ref => $count_file_ref,
      parameterFile2     => $basicMatrix_file,
      output_file_ext    => ".cpm.csv",
      use_tmp_folder     => 0,
      sh_direct          => 1,
      pbs                => {
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    };
    push( @$tasks, $cpm_task );

    my $deconvolution_task = "deconvolution_2_call";
    $config->{$deconvolution_task} = {
      class        => "CQS::ProgramWrapperOneToOne",
      perform      => 1,
      target_dir   => "$target_dir/$deconvolution_task",
      init_command => "
source activate " . getValue( $def, "deconvolution_conda" ) . "

export NUMEXPR_MAX_THREADS=12

",
      option                => "-b $basicMatrix_file -n __NAME__",
      interpretor           => "python3",
      check_program         => 1,
      program               => "../Deconvolution/do_deconvolve.py",
      source_ref            => $source_ref,
      source_arg            => "",
      parameterFile1_arg    => "-i",
      parameterFile1_ref    => [$cpm_task],
      output_to_same_folder => 1,
      no_input              => 1,
      no_output             => 1,
      output_arg            => "",
      output_file_ext       => "__NAME__/__NAME___deconvolutionCoefs.csv,__NAME__/__NAME__-NUSVR_supportVectors.csv",
      no_docker             => 1,
      use_tmp_folder        => 0,
      sh_direct             => 0,
      pbs                   => {
        "nodes"    => "1:ppn=12",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    };
    push( @$tasks, $deconvolution_task );

    my $merge_task = "deconvolution_3_merge";
    $config->{$merge_task} = {
      class                    => "CQS::ProgramWrapper",
      perform                  => 1,
      target_dir               => "$target_dir/$merge_task",
      init_command             => "source activate " . getValue( $def, "deconvolution_conda" ),
      option                   => "-i __FILE__ -o __NAME__",
      interpretor              => "python3",
      check_program            => 1,
      program                  => "../Deconvolution/merge.py",
      parameterSampleFile1_ref => $deconvolution_task,
      output_arg               => "-o",
      output_ext               => "",
      output_file_ext          => ".fractions.csv,.bestCoef.csv",
      no_docker                => 1,
      use_tmp_folder           => 0,
      sh_direct                => 0,
      pbs                      => {
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    };

    push( @$tasks, $merge_task );
  } ## end if ( $def->{perform_deconvolution...})

  my $deseq2taskname;
  my $webgestaltTaskName;
  my $gseaTaskName;
  my $linkTaskName;
  my $webgestaltHeatmapTaskName;
  if ( defined $def->{pairs} ) {
    my $de_prefix;
    my $de_source = $count_file_ref;
    if ( $def->{perform_proteincoding_gene} ) {
      if ( $def->{perform_combatseq} ) {
        $de_prefix = "combatseq_proteincoding_genetable";
        if ( defined $config->{genetable} ) {
          $de_source = [ "genetable_combatseq", ".proteincoding.count\$" ];
        }
      } ## end if ( $def->{perform_combatseq...})
      else {
        $de_prefix = "proteincoding_genetable";
        if ( defined $config->{genetable} ) {
          $de_source = [ "genetable", ".proteincoding.count\$" ];
        }
      } ## end else [ if ( $def->{perform_combatseq...})]
    } ## end if ( $def->{perform_proteincoding_gene...})
    else {
      if ( $def->{perform_combatseq} ) {
        $de_prefix = "combatseq_genetable";
      }
      else {
        $de_prefix = "genetable";
      }
    } ## end else [ if ( $def->{perform_proteincoding_gene...})]
    $deseq2taskname = addDEseq2( $config, $def, $tasks, $de_prefix, $de_source, $def->{target_dir}, $def->{DE_min_median_read} );

    if ( getValue( $def, "perform_webgestalt" ) ) {
      $webgestaltTaskName = addWebgestalt( $config, $def, $tasks, $target_dir, $deseq2taskname, [ $deseq2taskname, "sig_genename.txt\$" ] );

      $linkTaskName = $webgestaltTaskName . "_link_deseq2";
      if ( getValue( $def, "perform_link_webgestalt_deseq2_v2" ) ) {
        $config->{$linkTaskName} = {
          class                    => "CQS::UniqueRmd",
          perform                  => 1,
          target_dir               => $target_dir . "/" . getNextFolderIndex($def) . $linkTaskName,
          report_rmd_file          => "../Annotation/WebGestaltDeseq2.v2.rmd",
          additional_rmd_files     => "../CQS/reportFunctions.R;../Annotation/WebGestaltReportFunctions.r;../Annotation/WebGestaltDeseq2.v2.sub.rmd",
          option                   => "",
          parameterSampleFile1_ref => [ $webgestaltTaskName, ".txt\$" ],
          parameterSampleFile2_ref => [ $deseq2taskname,     "sig.csv\$" ],
          parameterSampleFile3     => {
            task_name   => getValue( $def, "task_name" ),
            email       => getValue( $def, "email" ),
            affiliation => $def->{"affiliation"},
          },
          output_file_ext          => ".webgestalt.html",
          output_other_ext         => ".webgestalt.rds",
          can_result_be_empty_file => 0,
          sh_direct                => 1,
          pbs                      => {
            "nodes"    => "1:ppn=1",
            "walltime" => "2",
            "mem"      => "10gb"
          },
        };
      } ## end if ( getValue( $def, "perform_link_webgestalt_deseq2_v2"...))
      else {
        $config->{$linkTaskName} = {
          class                       => "CQS::UniqueR",
          perform                     => 1,
          target_dir                  => $target_dir . "/" . getNextFolderIndex($def) . $linkTaskName,
          rtemplate                   => "../Annotation/WebGestaltReportFunctions.r;../Annotation/WebGestaltDeseq2.r",
          rReportTemplate             => "../Annotation/WebGestaltDeseq2.rmd",
          output_to_result_directory  => 1,
          output_perSample_file       => "parameterSampleFile1",
          output_perSample_file_regex => "enrichment_results_(.+).txt",
          output_perSample_file_ext   => ".html;.html.rds",
          parameterSampleFile1_ref    => [ $webgestaltTaskName, ".txt\$" ],
          parameterSampleFile2_ref    => [ $deseq2taskname,     "sig.csv\$" ],
          sh_direct                   => 1,
          rCode                       => "",
          pbs                         => {
            "nodes"    => "1:ppn=1",
            "walltime" => "23",
            "mem"      => "10gb"
          },
        };
      } ## end else [ if ( getValue( $def, "perform_link_webgestalt_deseq2_v2"...))]
      #if ( defined $def->{perform_link_webgestalt_deseq2} ) {
      push( @$tasks, $linkTaskName );

      if ( getValue( $def, "perform_webgestaltHeatmap" ) ) {
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
          parameterSampleFile2_ref   => [ $deseq2taskname,     "sig.csv\$" ],
          parameterSampleFile3_ref   => [ $deseq2taskname,     "vsd.csv\$" ],
          parameterSampleFile4_ref   => [ $deseq2taskname,     ".design\$" ],
          sh_direct                  => 1,
          rCode                      => "",
          pbs                        => {
            "nodes"    => "1:ppn=1",
            "walltime" => "23",
            "mem"      => "10gb"
          },
        };
        push( @$tasks, $webgestaltHeatmapTaskName );
      } ## end if ( getValue( $def, "perform_webgestaltHeatmap"...))
      #}
    } ## end if ( getValue( $def, "perform_webgestalt"...))

    #print(Dumper($def));

    if ( getValue( $def, "perform_gsea" ) ) {
      $gseaTaskName = $deseq2taskname . "_GSEA";

      if ( getValue( $def, "use_mouse_gsea_db", 0 ) ) {
        $gseaTaskName = $gseaTaskName . "_Mm";
      }
      else {
        $gseaTaskName = $gseaTaskName . "_Hs";
      }

      my $pairs = $config->{pairs};
      my $keys  = [ keys %$pairs ];
      #my $suffix = getDeseq2Suffix($config, $def, $deseq2taskname);

      add_gsea( $config, $def, $tasks, $target_dir, $gseaTaskName, [ $deseq2taskname, "_GSEA.rnk\$" ], $keys, "" );
    } ## end if ( getValue( $def, "perform_gsea"...))

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
          "nodes"    => "1:ppn=1",
          "walltime" => "23",
          "mem"      => "10gb"
        },
      };
      push( @$tasks, "keggprofile" );
    } ## end if ( $def->{perform_keggprofile...})
  } ## end if ( defined $def->{pairs...})

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
        "nodes"    => "1:ppn=" . $def->{max_thread},
        "walltime" => "23",
        "mem"      => "40gb"
      },
    };
    push( @$tasks, "rnaseqc" );
  } ## end if ( $def->{perform_rnaseqc...})

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
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "40gb"
      },
    };
    push( @$tasks, $rnaseqBamQC_task );

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
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      },
    };
    push( @$tasks, $rnaseqBamQCsummary_task );
  } ## end if ( $def->{perform_rnaseqBamQC...})

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
        "nodes"    => "1:ppn=8",
        "walltime" => "23",
        "mem"      => "40gb"
      },
    };
    push( @$tasks, "qc3" );
  } ## end if ( $def->{perform_qc3bam...})

  if ( $def->{perform_bamplot} ) {
    defined $def->{dataset_name} or die "Define dataset_name for bamplot first!";
    if ( not defined $def->{bamplot_gff} ) {
      defined $def->{gene_names} or die "Define gene_names for bamplot first, seperate by blank space!";
      defined $def->{add_chr}    or die "Define add_chr for bamplot first, check your genome sequence!";
    }
    add_bamplot( $config, $def, $tasks, $target_dir, $source_ref );
  } ## end if ( $def->{perform_bamplot...})

  if ( $def->{perform_call_variants} ) {
    my $fasta = getValue( $def, "fasta_file" );
    my $dbsnp = getValue( $def, "dbsnp" );

    my $gatk_index     = $def;
    my $gatk_index_snv = "SNV_Index";
    my $gatk_prefix;

    my $refine_task;
    my $pass_task;
    if ( getValue( $def, "perform_call_variants_by_gatk4", 1 ) ) {
      $refine_task = "gatk4_refine";
      $config->{$refine_task} = {
        class      => "CQS::ProgramWrapperOneToOne",
        target_dir => $target_dir . "/" . getNextFolderIndex($def) . $refine_task,
        option     => "

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
    gatk --java-options \"-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal -XX:+PrintGCDetails -Xloggc:gc_log.log -Xms4000m\" \\
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
      gatk --java-options \"-XX:+PrintFlagsFinal \\
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
        suffix                => "_rrf",
        docker_prefix         => "gatk4_",
        interpretor           => "",
        program               => "",
        check_program         => 0,
        source_arg            => "",
        source_ref            => $source_ref,
        output_arg            => "",
        output_file_ext       => ".rmdup.split.recal.bam",
        output_to_same_folder => 1,
        sh_direct             => getValue( $def, "sh_direct", 0 ),
        pbs                   => {
          "nodes"    => "1:ppn=1",
          "walltime" => "23",
          "mem"      => "40gb"
        }
      };

      push( @$tasks, $refine_task );

      $gatk_prefix = $refine_task . "_SNV_";

      my $hc_task = $gatk_prefix . getNextIndex( $gatk_index, $gatk_index_snv ) . "_hc";
      $config->{$hc_task} = {
        class         => "GATK4::HaplotypeCaller",
        perform       => 1,
        target_dir    => $target_dir . "/" . getNextFolderIndex($def) . $hc_task,
        option        => getValue( $def, "HaplotypeCaller_option", "--soft-clip-low-quality-ends true --dont-use-soft-clipped-bases true --standard-min-confidence-threshold-for-calling 20" ),
        source_ref    => $refine_task,
        java_option   => "",
        bed_file      => $def->{covered_bed},
        fasta_file    => $fasta,
        extension     => ".vcf.gz",
        by_chromosome => 0,                                                                                                                                                                       #since we have the bed file, we cannot use by_chromosome.
        gvcf          => 0,                                                                                                                                                                       #http://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq
        sh_direct     => 0,
        pbs           => {
          "nodes"    => "1:ppn=8",
          "walltime" => getValue( $def, "HaplotypeCaller_walltime", "23" ),
          "mem"      => "40gb"
        },
      };
      push( @$tasks, $hc_task );

      my $min_dp      = getValue( $def, "SNV_minimum_depth", 10 );
      my $filter_task = $gatk_prefix . getNextIndex( $gatk_index, $gatk_index_snv ) . "_filter";
      $config->{$filter_task} = {
        class      => "CQS::ProgramWrapperOneToOne",
        target_dir => $target_dir . "/" . $filter_task,
        option     => "
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
        suffix                       => "_vf",
        docker_prefix                => "gatk4_",
        interpretor                  => "",
        program                      => "",
        check_program                => 0,
        source_arg                   => "--V",
        source_ref                   => $hc_task,
        other_localization_ext_array => [".tbi"],
        output_arg                   => "-O",
        output_file_ext              => ".variant_filtered.vcf.gz",
        output_to_same_folder        => 1,
        sh_direct                    => getValue( $def, "sh_direct", 0 ),
        pbs                          => {
          "nodes"    => "1:ppn=1",
          "walltime" => "10",
          "mem"      => "10gb"
        }
      };
      push( @$tasks, $filter_task );

      my $lefttrim_task = $gatk_prefix . getNextIndex( $gatk_index, $gatk_index_snv ) . "_lefttrim";
      $config->{$lefttrim_task} = {
        class         => "GATK4::LeftTrim",
        perform       => 1,
        target_dir    => $target_dir . "/" . $lefttrim_task,
        source_ref    => $filter_task,
        option        => "",
        docker_prefix => "exome_",
        fasta_file    => $fasta,
        extension     => ".variant_filtered.norm.nospan.vcf.gz",
        sh_direct     => 0,
        pbs           => {
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      };
      push( @$tasks, $lefttrim_task );

      $pass_task = $gatk_prefix . getNextIndex( $gatk_index, $gatk_index_snv ) . "_merge";
      $config->{$pass_task} = {
        class      => "CQS::ProgramWrapper",
        target_dir => $target_dir . "/" . $pass_task,
        option     => "-i __FILE__ -o __NAME__.pass.vcf.gz

#__OUTPUT__
",
        suffix          => "_mg",
        interpretor     => "python3",
        docker_prefix   => "exome_",
        program         => "../GATK4/combineVCFs.py",
        check_program   => 1,
        source_arg      => "-i",
        source_ref      => $lefttrim_task,
        output_arg      => "-o",
        output_file_ext => ".pass.vcf.gz",
        sh_direct       => 1,
        pbs             => {
          "nodes"    => "1:ppn=1",
          "walltime" => "10",
          "mem"      => "10gb"
        }
      };
      push( @$tasks, $pass_task );
    } ## end if ( getValue( $def, "perform_call_variants_by_gatk4"...))
    else {
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
          "nodes"    => "1:ppn=8",
          "walltime" => "23",
          "mem"      => "40gb"
        },
      };
      push( @$tasks, $refine_task );

      my $hc_task = $refine_task . "_hc";
      $config->{$hc_task} = {
        class         => "GATK::HaplotypeCaller",
        perform       => 1,
        target_dir    => $target_dir . "/" . getNextFolderIndex($def) . $hc_task,
        option        => getValue( $def, "HaplotypeCaller_option", "--soft-clip-low-quality-ends true --dont-use-soft-clipped-bases true --standard-min-confidence-threshold-for-calling 20" ),
        source_ref    => $refine_task,
        java_option   => "",
        fasta_file    => $fasta,
        gatk_jar      => $gatk,
        extension     => ".vcf",
        by_chromosome => 0,                                                                                                                                                                       #since we have the bed file, we cannot use by_chromosome.
        gvcf          => 0,                                                                                                                                                                       #http://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq
        sh_direct     => 0,
        pbs           => {
          "nodes"    => "1:ppn=8",
          "walltime" => "23",
          "mem"      => "40gb"
        },
      };
      push( @$tasks, $hc_task );

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
          "nodes"    => "1:ppn=8",
          "walltime" => "23",
          "mem"      => "40gb"
        },
      };
      push( @$tasks, $filter_task );

      $pass_task = $filter_task;
    } ## end else [ if ( getValue( $def, "perform_call_variants_by_gatk4"...))]

    $multiqc_depedents = $refine_task;

    if ( $def->{filter_variants_by_allele_frequency} ) {
      my $maf_filter_task = $gatk_prefix . getNextIndex( $gatk_index, $gatk_index_snv ) . "_filterMAF";
      add_maf_filter( $config, $def, $tasks, $target_dir, $maf_filter_task, $pass_task );
      $pass_task = $maf_filter_task;
    }

    if ( $def->{perform_annovar} ) {
      my $annovar_task = addAnnovar( $config, $def, $tasks, $target_dir, $pass_task, undef, $gatk_prefix, $gatk_index, $gatk_index_snv );

      if ( $def->{annovar_param} =~ /exac/ || $def->{annovar_param} =~ /1000g/ || $def->{annovar_param} =~ /gnomad/ ) {
        my $annovar_filter_task = addAnnovarFilter( $config, $def, $tasks, $target_dir, $annovar_task, $gatk_prefix, $gatk_index, $gatk_index_snv );

        if ( defined $def->{annotation_genes} ) {
          addAnnovarFilterGeneannotation( $config, $def, $tasks, $target_dir, $annovar_filter_task );
        }

        addAnnovarMafReport( $config, $def, $tasks, $target_dir, $annovar_filter_task, $gatk_prefix, $gatk_index, $gatk_index_snv );
      } ## end if ( $def->{annovar_param...})
    } ## end if ( $def->{perform_annovar...})
  } ## end if ( $def->{perform_call_variants...})

  if ( getValue( $def, "perform_multiqc" ) ) {
    addMultiQC( $config, $def, $tasks, $target_dir, $target_dir, $multiqc_depedents );
  }

  if ( getValue( $def, "perform_report" ) ) {
    my @report_files = ();
    my @report_names = ();
    my @copy_files   = ();

    my $version_files = get_version_files($config);

    if ( defined $config->{fastqc_raw_summary} ) {
      push( @report_files, "fastqc_raw_summary", ".FastQC.baseQuality.tsv.png" );
      push( @report_files, "fastqc_raw_summary", ".FastQC.sequenceGC.tsv.png" );
      push( @report_files, "fastqc_raw_summary", ".FastQC.adapter.tsv.png" );
      push( @report_names, "fastqc_raw_per_base_sequence_quality", "fastqc_raw_per_sequence_gc_content", "fastqc_raw_adapter_content" );
    } ## end if ( defined $config->...)

    if ( defined $config->{fastqc_post_trim_summary} ) {
      push( @report_files, "fastqc_post_trim_summary", ".FastQC.baseQuality.tsv.png" );
      push( @report_files, "fastqc_post_trim_summary", ".FastQC.sequenceGC.tsv.png" );
      push( @report_files, "fastqc_post_trim_summary", ".FastQC.adapter.tsv.png" );
      push( @report_names, "fastqc_post_trim_per_base_sequence_quality", "fastqc_post_trim_per_sequence_gc_content", "fastqc_post_trim_adapter_content" );
    } ## end if ( defined $config->...)

    if ( defined $config->{fastqc_count_vis} ) {
      push( @report_files, "fastqc_count_vis", ".countInFastQcVis.Result.png" );
      push( @report_names, "fastqc_post_trim_reads" );
    }

    if ( defined $config->{star_featurecount_summary} ) {
      push( @report_files, "star_featurecount_summary", ".STARSummary.csv.png" );
      push( @report_files, "star_featurecount_summary", ".STARSummary.csv\$" );
      push( @report_names, "STAR_summary",              "STAR_summary_table" );

      push( @report_files, "star_featurecount_summary", ".chromosome.png\$" );
      push( @report_files, "star_featurecount_summary", ".chromosome.csv\$" );
      push( @report_names, "STAR_chromosome_png",       "STAR_chromosome_table" );

      push( @report_files, "star_featurecount_summary", ".FeatureCountSummary.csv.png\$" );
      push( @report_files, "star_featurecount_summary", ".FeatureCountSummary.csv\$" );
      push( @report_names, "featureCounts_table_png",   "featureCounts_table" );

      push( @report_files, "star_featurecount_summary", ".gene.count.png\$" );
      push( @report_files, "star_featurecount_summary", ".gene.count.csv\$" );
      push( @report_names, "featureCounts_gene_png",    "featureCounts_gene_table" );
    } ## end if ( defined $config->...)

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
      if ( $def->{perform_proteincoding_gene} ) {
        push( @report_files, "genetable", ".fpkm.proteincoding.tsv\$" );
      }
      else {
        push( @report_files, "genetable", ".fpkm.tsv\$" );
      }
      push( @report_names, "genetable_fpkm" );
    } ## end if ( defined $config->...)

    if ( defined $config->{genetable_combatseq} ) {
      push( @copy_files, "genetable_combatseq", ".count\$" );
    }

    if ( defined $config->{star_fusion_summary} ) {
      push( @report_files, "star_fusion_summary", ".fusion_count.png" );
      push( @report_names, "star_fusion_png" );
    }

    my @corr_tasks = ( "genetable_correlation", "genetable_combatseq_correlation" );
    foreach my $cur_corr_task (@corr_tasks) {
      if ( defined $config->{$cur_corr_task} ) {
        my $suffix = $config->{$cur_corr_task}{suffix};
        if ( !defined $suffix ) {
          $suffix = "";
        }
        my $pcoding = $def->{perform_proteincoding_gene} ? ".proteincoding.count" : ".count";

        my $titles = { "all" => "" };
        if ( is_not_array($count_file_ref) ) {    #count file directly
          $titles->{all} = basename($count_file_ref);
          $pcoding = "";
        }
        if ( defined $config->{$cur_corr_task}{parameterSampleFile2} ) {
          my $correlationGroups = $config->{$cur_corr_task}{parameterSampleFile2};
          for my $correlationTitle ( keys %$correlationGroups ) {
            my $groups = $correlationGroups->{$correlationTitle};
            if ( is_hash($groups) ) {
              if ( $correlationTitle ne "all" ) {
                $correlationTitle =~ s/\\s+/_/g;
                $titles->{$correlationTitle} = "." . $correlationTitle;
              }
            } ## end if ( is_hash($groups) )
          } ## end for my $correlationTitle...
        } ## end if ( defined $config->...)

        for my $title ( keys %$titles ) {
          my $title_prefix = defined $config->{genetable_combatseq} ? $cur_corr_task . "." . $title : $title;
          push( @report_files, $cur_corr_task, $pcoding . $suffix . $titles->{$title} . ".density.png", $cur_corr_task, $pcoding . $suffix . $titles->{$title} . ".heatmap.png", $cur_corr_task, $pcoding . $suffix . $titles->{$title} . ".PCA.png" );
          push( @report_names, $title_prefix . "_correlation_density", $title_prefix . "_correlation_heatmap", $title_prefix . "_correlation_PCA" );
        }
      } ## end if ( defined $config->...)
    } ## end foreach my $cur_corr_task (...)

    if ( ( defined $deseq2taskname ) && ( defined $config->{$deseq2taskname} ) ) {
      my $suffix = getDeseq2Suffix( $config, $def, $deseq2taskname );

      my $pairs = $config->{pairs};

      if ( scalar( keys %$pairs ) > 1 ) {
        push( @report_files, $deseq2taskname, "/" . $taskName . ".define.*DESeq2_volcanoPlot.png" );
        push( @report_names, "deseq2_volcano_plot" );
      }
      else {
        if ( $def->{REPORT_use_enhanced_volcano} ) {
          push( @report_files, $deseq2taskname, "_DESeq2_volcanoEnhanced.png" );
        }
        else {
          push( @report_files, $deseq2taskname, "_DESeq2_volcanoPlot.png" );
        }
        push( @report_names, "deseq2_volcano_plot" );
      } ## end else [ if ( scalar( keys %$pairs...))]
      for my $key ( keys %$pairs ) {
        push( @report_files, $deseq2taskname, "/" . $key . $suffix . "_DESeq2_sig.csv" );
        push( @report_names, "deseq2_" . $key );

        push( @report_files, $deseq2taskname, "/" . $key . ".design" );
        push( @report_names, "deseq2_" . $key . "_design" );

        push( @report_files, $deseq2taskname, "/" . $key . $suffix . "_geneAll_DESeq2-vsd-heatmap.png" );
        push( @report_names, "deseq2_" . $key . "_heatmap" );

        push( @report_files, $deseq2taskname, "/" . $key . $suffix . "_geneAll_DESeq2-vsd-pca.png" );
        push( @report_names, "deseq2_" . $key . "_pca" );
      } ## end for my $key ( keys %$pairs)
      push( @copy_files, $deseq2taskname, "_DESeq2.csv" );
      push( @copy_files, $deseq2taskname, "_DESeq2_sig.csv" );
      push( @copy_files, $deseq2taskname, "_DESeq2-vsd.csv" );

      #push( @copy_files, $deseq2taskname, "_DESeq2_GSEA.rnk" );
      #push( @copy_files, $deseq2taskname, "_DESeq2_sig_genename.txt" );
      #push( @copy_files, $deseq2taskname, "heatmap.png" );
      #push( @copy_files, $deseq2taskname, "pca.pdf" );
    } ## end if ( ( defined $deseq2taskname...))

    my $hasFunctionalEnrichment = 0;
    if ( defined $webgestaltTaskName ) {
      push( @copy_files, $webgestaltTaskName, "_geneontology_Biological_Process\$" );
      push( @copy_files, $webgestaltTaskName, "_geneontology_Cellular_Component\$" );
      push( @copy_files, $webgestaltTaskName, "_geneontology_Molecular_Function\$" );
      push( @copy_files, $webgestaltTaskName, "_pathway_KEGG\$" );

      if ( defined $linkTaskName && defined $config->{$linkTaskName} ) {
        push( @copy_files, $linkTaskName, ".html\$" );

        push( @report_files, $linkTaskName, ".rds" );
        push( @report_names, "WebGestalt_deseq2" );
      } ## end if ( defined $linkTaskName...)
      else {
        my $pairs = $config->{pairs};
        for my $key ( keys %$pairs ) {
          push( @report_files, $webgestaltTaskName, "enrichment_results_" . $key . "_geneontology_Biological_Process.txt" );
          push( @report_files, $webgestaltTaskName, "enrichment_results_" . $key . "_geneontology_Cellular_Component.txt" );
          push( @report_files, $webgestaltTaskName, "enrichment_results_" . $key . "_geneontology_Molecular_Function.txt" );
          push( @report_files, $webgestaltTaskName, "enrichment_results_" . $key . "_pathway_KEGG.txt" );

          push( @report_names, "WebGestalt_GO_BP_" . $key );
          push( @report_names, "WebGestalt_GO_CC_" . $key );
          push( @report_names, "WebGestalt_GO_MF_" . $key );
          push( @report_names, "WebGestalt_KEGG_" . $key );
        } ## end for my $key ( keys %$pairs)
      } ## end else [ if ( defined $linkTaskName...)]
      $hasFunctionalEnrichment = 1;
    } ## end if ( defined $webgestaltTaskName)

    if ( defined $gseaTaskName ) {
      push( @copy_files, $gseaTaskName, ".gsea\$" );
      #my $suffix = getDeseq2Suffix($config, $def, $deseq2taskname);

      my $pairs = $config->{pairs};
      for my $key ( keys %$pairs ) {
        #push( @report_files, $gseaTaskName, "/" . $key . $suffix . "_.*gsea.csv" );
        push( @report_files, $gseaTaskName, "/" . $key . ".gsea.csv" );
        push( @report_names, "gsea_" . $key );
      }

      my $gsea_report_task = $gseaTaskName . "_report";
      push( @report_files, $gsea_report_task, "gsea_files.csv" );
      push( @report_names, "report_gsea" );

      $hasFunctionalEnrichment = 1;
    } ## end if ( defined $gseaTaskName)

    my $fcOptions      = getValue( $def, "featureCount_option" );
    my $fcMultiMapping = ( $fcOptions =~ /-m/ ) ? "TRUE" : "FALSE";
    my $options        = {
      "DE_min_median_read"                 => [ getValue( $def, "DE_min_median_read" ) ],
      "perform_proteincoding_gene"         => [ getValue( $def, "perform_proteincoding_gene" ) ? "TRUE" : "FALSE" ],
      "DE_fold_change"                     => [ getValue( $def, "DE_fold_change",    2 ) ],
      "DE_pvalue"                          => [ getValue( $def, "DE_pvalue",         0.05 ) ],
      "DE_use_raw_pvalue"                  => [ getValue( $def, "DE_use_raw_pvalue", 0 ) ],
      "DE_combatseq"                       => [ getValue( $def, "DE_combatseq",      0 ) ],
      "featureCounts_UseMultiMappingReads" => [$fcMultiMapping],
      "top25cv_in_hca"                     => [ getValue( $def, "top25cv_in_hca" ) ? "TRUE" : "FALSE" ],
      "task_name"                          => $taskName,
      "out.width"                          => getValue( $def, "report.out.width",    "80%" ),
      "remove_chrM_genes"                  => getValue( $def, "remove_chrM_genes",   0 ),
      "adapter"                            => getValue( $def, "adapter",             "" ),
      "cutadapt_option"                    => getValue( $def, "cutadapt_option",     "" ),
      "featureCount_option"                => getValue( $def, "featureCount_option", "" )
    };

    if ( $def->{introduction_rmd} ) {
      $options->{introduction_rmd} = $def->{introduction_rmd};
    }

    $config->{report} = {
      class                      => "CQS::BuildReport",
      perform                    => 1,
      target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . "report",
      report_rmd_file            => "../Pipeline/RNASeq.Rmd",
      additional_rmd_files       => "../Pipeline/Pipeline.R;reportFunctions.R;../Annotation/WebGestaltDeseq2.v2.sub.rmd",
      parameterSampleFile1_ref   => \@report_files,
      parameterSampleFile1_names => \@report_names,
      parameterSampleFile2       => $options,
      parameterSampleFile3_ref   => \@copy_files,
      parameterSampleFile4       => $version_files,
      parameterSampleFile5       => $def->{software_version},
      parameterSampleFile6       => $def->{groups},

      sh_direct => 1,
      pbs       => {
        "nodes"    => "1:ppn=1",
        "walltime" => "1",
        "mem"      => "20gb"
      },
    };
    push( @$tasks, "report" );
  } ## end if ( getValue( $def, "perform_report"...))

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => { tasks => $tasks, },
    sh_direct  => 0,
    cluster    => $cluster,
    pbs        => {
      "nodes"    => "1:ppn=" . $def->{max_thread},
      "walltime" => $def->{sequencetask_run_time},
      "mem"      => "40gb"
    },
  };

  return ($config);
} ## end sub getRNASeqConfig


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
} ## end sub performRNASeq


sub performRNASeqTask {
  my ( $def, $task ) = @_;

  my $config = getRNASeqConfig($def);

  performTask( $config, $task );

  return $config;
} ## end sub performRNASeqTask

1;
