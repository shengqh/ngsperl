#!/usr/bin/perl
package Pipeline::RNASeq;

use strict;
use warnings;
use List::Util qw(first);
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

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "ALL" );
  initDefaultValue( $def, "cluster", "slurm");

  initDefaultValue( $def, "perform_preprocessing", 1 );
  initDefaultValue( $def, "perform_mapping",       1 );
  initDefaultValue( $def, "perform_counting",      1 );
  initDefaultValue( $def, "perform_count_table",   1 );
  initDefaultValue( $def, "perform_correlation",   1 );
  initDefaultValue( $def, "perform_rnaseqc",       0 );
  initDefaultValue( $def, "perform_qc3bam",        0 );
  initDefaultValue( $def, "perform_bamplot",       0 );
  initDefaultValue( $def, "perform_call_variants", 0 );
  initDefaultValue( $def, "perform_multiqc",       1 );
  initDefaultValue( $def, "perform_webgestalt",    0 );
  initDefaultValue( $def, "perform_gsea",          0 );
  initDefaultValue( $def, "perform_report",        1 );

  initDefaultValue( $def, "perform_cutadapt", 0 );

  initDefaultValue( $def, "featureCount_option",        "-g gene_id -t exon" );
  initDefaultValue( $def, "aligner",                    "star" );
  initDefaultValue( $def, "star_option",                "--twopassMode Basic --outSAMprimaryFlag AllBestScore" );
  initDefaultValue( $def, "use_pearson_in_hca",         1 );
  initDefaultValue( $def, "top25cv_in_hca",             0 );
  initDefaultValue( $def, "use_green_red_color_in_hca", 1 );
  initDefaultValue( $def, "output_bam_to_same_folder",  1 );
  initDefaultValue( $def, "show_label_PCA",  1 );
  
  initDefaultValue( $def, "max_thread",                 8 );
  initDefaultValue( $def, "sequencetask_run_time",      '24' );

  initDefaultValue( $def, "perform_keggprofile",      0 );
  initDefaultValue( $def, "keggprofile_useRawPValue", 0 );
  initDefaultValue( $def, "pairend",                  1 );

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
  initDefaultValue( $def, "perform_DE_proteincoding_gene",   0 );
  initDefaultValue( $def, "perform_proteincoding_gene",      getValue($def, "perform_DE_proteincoding_gene") );
  
  initDefaultValue( $def, "outputPdf",      0 );
  initDefaultValue( $def, "outputPng",      1 );
  
  return $def;
}

sub getRNASeqConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

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

  my $target_dir      = $def->{target_dir};
  my $groups_ref      = defined $def->{groups} ? "groups" : undef;
  my $aligner         = $def->{aligner};
  my $star_option     = $def->{star_option};
  my $count_table_ref = "files";
  
  my $multiqc_depedents = $source_ref;
  
  my $count_file_ref = $def->{count_file};
  if ( $def->{perform_mapping} && $def->{perform_counting} && ( $aligner eq "star" ) && $def->{perform_star_featurecount} ) {
    my $aligner_index   = $def->{star_index} or die "Define star_index at definition first";
    my $starFolder      = $target_dir . "/" . getNextFolderIndex($def) . "star_featurecount";
    my $transcript_gtf  = $def->{transcript_gtf} or die "Define transcript_gtf at definition first";
    my $configAlignment = {
      "star_featurecount" => {
        class                     => "Alignment::STARFeatureCount",
        perform                   => 1,
        target_dir                => $starFolder,
        option                    => $star_option,
        source_ref                => $source_ref,
        genome_dir                => $aligner_index,
        output_sort_by_coordinate => 1,
        output_to_same_folder     => $def->{output_bam_to_same_folder},
        featureCount_option       => getValue( $def, "featureCount_option" ),
        star_location             => $def->{star_location},
        gff_file                  => $transcript_gtf,
        ispairend                 => getValue( $def, "pairend" ),
        sh_direct                 => 0,
        pbs                       => {
          "email"     => $email,
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=" . $def->{max_thread},
          "walltime"  => "72",
          "mem"       => "40gb"
        },
      },
      "star_summary" => {
        class      => "Alignment::STARSummary",
        perform    => 1,
        target_dir => $starFolder,
        option     => "",
        source_ref => [ "star_featurecount", "_Log.final.out" ],
        sh_direct  => 1,
        pbs        => {
          "email"     => $email,
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "72",
          "mem"       => "40gb"
        },
      },
    };

    $source_ref      = [ "star_featurecount", "_Aligned.sortedByCoord.out.bam\$" ];
    $count_table_ref = [ "star_featurecount", ".count\$" ];
    push @$individual, ("star_featurecount");
    push @$summary,    ("star_summary");
    $config = merge( $config, $configAlignment );
    
    $multiqc_depedents = "star_featurecount";
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
              "email"     => $email,
              "emailType" => $def->{emailType},
              "nodes"     => "1:ppn=" . $def->{max_thread},
              "walltime"  => "24",
              "mem"       => "40gb"
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
              "email"     => $email,
              "emailType" => $def->{emailType},
              "nodes"     => "1:ppn=1",
              "walltime"  => "2",
              "mem"       => "10gb"
            },
          }
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
              "email"     => $email,
              "emailType" => $def->{emailType},
              "nodes"     => "1:ppn=" . $def->{max_thread},
              "walltime"  => "24",
              "mem"       => "40gb"
            },
          },
        };
        $source_ref = [ "hisat2", ".bam\$" ];
      }

      $config = merge( $config, $configAlignment );
      push @$individual, $aligner;
      $multiqc_depedents = "hisat2";
    }

    if ( $def->{perform_counting} ) {
      my $cqstools       = $def->{cqstools};
      my $transcript_gtf = $def->{transcript_gtf} or die "Define transcript_gtf at definition first";
      if ( $def->{additional_bam_files} ) {
        push @$source_ref, "additional_bam_files";
      }
      $config->{"featurecount"} = {
        class      => "Count::FeatureCounts",
        perform    => 1,
        target_dir => $target_dir . "/" . getNextFolderIndex($def) . "featurecount",
        option     => "-g gene_id -t exon",
        source_ref => $source_ref,
        gff_file   => $transcript_gtf,
        ispairend  => 1,
        sh_direct  => 0,
        pbs        => {
          "email"     => $email,
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "24",
          "mem"       => "40gb"
        },
      };

      push @$individual, "featurecount";
      $count_table_ref = [ "featurecount", ".count\$" ];
      $multiqc_depedents = "featurecount";
    }
  }

  my $perform_count_table = $def->{perform_counting} || $def->{perform_count_table};

  if ($perform_count_table) {
    my $cqstools      = $def->{cqstools};
    my $name_map_file = $def->{name_map_file};
    $config->{"genetable"} = {
      class                     => "CQS::CQSDatatable",
      perform                   => 1,
      target_dir                => $target_dir . "/" . getNextFolderIndex($def) . "genetable",
      option                    => "-k 0 -v 6 -e --fillMissingWithZero",
      source_ref                => $count_table_ref,
      output_proteincoding_gene => $def->{perform_proteincoding_gene},
      name_map_file             => $name_map_file,
      cqs_tools                 => $cqstools,
      sh_direct                 => 1,
      pbs                       => {
        "email"     => $email,
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "24",
        "mem"       => "10gb"
      },
    };

    push @$summary, "genetable";

    $count_file_ref = [ "genetable", "(?<!proteincoding).count\$" ];
    if ( $def->{perform_proteincoding_gene} ) {
      push @$count_file_ref, "genetable", ".proteincoding.count\$";
    }
  }

  if ( $def->{perform_correlation} ) {
    my $cor_dir   = ( defined $config->{genetable} ) ? $config->{genetable}{target_dir} : $target_dir . "/" . getNextFolderIndex($def) . "genetable_correlation";
    my $gene_file = $def->{correlation_gene_file};
    my $rCode     = getOutputFormat($def);
    if ( defined $gene_file ) {
      $rCode = $rCode . "suffix<-\"_genes\"; ";
    }
    
    $config->{"genetable_correlation"} = {
      class       => "CQS::UniqueR",
      perform     => 1,
      suffix      => "_cor",
      rCode       => $rCode . "usePearsonInHCA<-" . $def->{use_pearson_in_hca} . "; useGreenRedColorInHCA<-" . $def->{use_green_red_color_in_hca} . "; top25cvInHCA<-" . $def->{top25cv_in_hca} . "; ",
      target_dir  => $cor_dir,
      rtemplate   => "countTableVisFunctions.R,countTableGroupCorrelation.R",
      output_file => "parameterSampleFile1",
      output_file_ext          => ".Correlation.png;.density.png;.heatmap.png;.PCA.png;.Correlation.Cluster.png",
      parameterFile1           => $gene_file,
      parameterSampleFile2_ref => $groups_ref,
      sh_direct                => 1,
      pbs                      => {
        "email"     => $email,
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "24",
        "mem"       => "20gb"
      },
    };
    if ( ref($count_file_ref) eq "ARRAY" ) {
      $config->{genetable_correlation}{parameterSampleFile1_ref} = $count_file_ref;
    }
    else {
      $config->{genetable_correlation}{parameterSampleFile1} = $count_file_ref;
    }
    push @$summary, "genetable_correlation";
  }

  my $deseq2taskname;
  my $webgestaltTaskName;
  my $gseaTaskName;
  if ( defined $def->{pairs} ) {
    if ( $def->{perform_proteincoding_gene} ) {
      $deseq2taskname = addDEseq2( $config, $def, $summary, "proteincoding_genetable", [ "genetable", ".proteincoding.count\$" ], $def->{target_dir}, $def->{DE_min_median_read} );
    }
    else {
      $deseq2taskname = addDEseq2( $config, $def, $summary, "genetable", $count_file_ref, $def->{target_dir}, $def->{DE_min_median_read} );
    }

    if ( getValue( $def, "perform_webgestalt" ) ) {
      $webgestaltTaskName = $deseq2taskname . "_WebGestalt";
      $config->{$webgestaltTaskName} = {
        class            => "Annotation::WebGestaltR",
        perform          => 1,
        target_dir       => $target_dir . "/" . getNextFolderIndex($def) . $webgestaltTaskName,
        option           => "",
        source_ref       => [ $deseq2taskname, "sig_genename.txt\$" ],
        organism         => getValue( $def, "webgestalt_organism" ),
        interestGeneType => $def->{interestGeneType},
        referenceSet     => $def->{referenceSet},
        sh_direct        => 1,
        pbs              => {
          "email"     => $email,
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "24",
          "mem"       => "10gb"
        },
      };
      push @$summary, "$webgestaltTaskName";
    }

    if ( getValue( $def, "perform_gsea" ) ) {
      my $gsea_jar        = $def->{gsea_jar}        or die "Define gsea_jar at definition first";
      my $gsea_db         = $def->{gsea_db}         or die "Define gsea_db at definition first";
      my $gsea_categories = $def->{gsea_categories} or die "Define gsea_categories at definition first";

      $gseaTaskName = $deseq2taskname . "_GSEA";

      #my $gseaCategories = "'h.all.v6.1.symbols.gmt','c2.all.v6.1.symbols.gmt','c5.all.v6.1.symbols.gmt','c6.all.v6.1.symbols.gmt','c7.all.v6.1.symbols.gmt'";
      $config->{$gseaTaskName} = {
        class                      => "CQS::UniqueR",
        perform                    => 1,
        target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . $gseaTaskName,
        rtemplate                  => "GSEAPerform.R",
        rReportTemplate            => "GSEAReport.Rmd",
        output_to_result_directory => 1,
        output_file                => "parameterSampleFile1",
        output_file_ext            => ".gsea.html;.gsea.csv;.gsea;",
        parameterSampleFile1_ref   => [ $deseq2taskname, "_GSEA.rnk\$" ],
        sh_direct                  => 1,
        rCode                      => "gseaDb='" . $gsea_db . "'; gseaJar='" . $gsea_jar . "'; gseaCategories=c(" . $gsea_categories . "); makeReport=F;",
        pbs                        => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "24",
          "mem"       => "10gb"
        },
      };
      push( @$summary, $gseaTaskName );
    }

    if ( $def->{perform_keggprofile} ) {
      my $keggprofile_useRawPValue = $def->{keggprofile_useRawPValue} or die "Define keggprofile_useRawPValue at definition first";
      $config->{keggprofile} = {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $target_dir . "/" . getNextFolderIndex($def) . "keggprofile",
        rtemplate                => "KEGGprofilePerform.R",
        output_file              => "",
        output_file_ext          => ".KEGG.csv",
        parameterSampleFile1_ref => [ $deseq2taskname, "_DESeq2.csv\$" ],
        sh_direct                => 1,
        rCode                    => "useRawPValue='" . $keggprofile_useRawPValue . "';",
        pbs                      => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "24",
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
        "email"     => $email,
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=" . $def->{max_thread},
        "walltime"  => "72",
        "mem"       => "40gb"
      },
    };
    push( @$summary, "rnaseqc" );
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
        "email"     => $email,
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "72",
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
          "email"     => $email,
          "emailType" => $def->{emailType},
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
          "email"     => $email,
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "24",
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
          "email"     => $email,
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "24",
          "mem"       => "10gb"
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
        "email"     => $email,
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=8",
        "walltime"  => "72",
        "mem"       => "40gb"
      },
    };
    
    $multiqc_depedents = "refine";

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
        "email"     => $email,
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=8",
        "walltime"  => "72",
        "mem"       => "40gb"
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
        "email"     => $email,
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=8",
        "walltime"  => "72",
        "mem"       => "40gb"
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
        "email"     => $email,
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "24",
        "mem"       => "10gb"
      },
    };
    push( @$individual, "refine",           "refine_hc" );
    push( @$summary,    "refine_hc_filter", "refine_hc_filter_annovar" );

  }

  if ( getValue( $def, "perform_multiqc" ) ) {
    addMultiQC( $config, $def, $summary, $target_dir, $target_dir, $multiqc_depedents );
  }

  if ( getValue( $def, "perform_report" ) ) {
    my @report_files = ();
    my @report_names = ();
    my @copy_files   = ();
    if ( defined $config->{multiqc} ) {
      my @configKeys = keys %$config;
      if ( first { $_ =~ 'fastqc' } @configKeys ) {
        push( @report_files, "multiqc",                          "fastqc_per_base_sequence_quality_plot_1.png" );
        push( @report_files, "multiqc",                          "fastqc_per_sequence_gc_content_plot_Percentages.png" );
        push( @report_files, "multiqc",                          "fastqc_adapter_content_plot_1.png" );
        push( @report_names, "fastqc_per_base_sequence_quality", "fastqc_per_sequence_gc_content", "fastqc_adapter_content" );
      }

      if ( defined $config->{star_featurecount} || defined $config->{featurecount} ) {
        push( @report_files, "multiqc", "multiqc_featureCounts.txt" );
        push( @report_names, "featureCounts_table" );
      }
    }

    if ( defined $config->{star_summary} ) {
      push( @report_files, "star_summary", ".STARSummary.csv.png" );
      push( @report_files, "star_summary", ".STARSummary.csv\$" );
      push( @report_names, "STAR_summary", "STAR_summary_table" );
    }

    push( @copy_files, "genetable", ".count\$", "genetable", ".fpkm.tsv" );

    if ( defined $config->{genetable_correlation} ) {
      my $pcoding = $def->{perform_proteincoding_gene} ? ".proteincoding.count" : "";
      push( @report_files,
        "genetable_correlation", $pcoding . ".density.png",
        "genetable_correlation", $pcoding . ".heatmap.png",
        "genetable_correlation", $pcoding . ".PCA.png",
        "genetable_correlation", $pcoding . ".Correlation.Cluster.png" );
      push( @report_names, "correlation_density", "correlation_heatmap", "correlation_PCA", "correlation_cluster" );
    }

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

      my $pairs = $config->{pairs};

      if ( scalar( keys %$pairs ) > 1 ) {
        push( @report_files, $deseq2taskname, $taskName . ".define.*DESeq2_volcanoPlot.png" );
        push( @report_names, "deseq2_volcano_plot" );
      }
      else {
        push( @report_files, $deseq2taskname, "_DESeq2_volcanoPlot.png" );
        push( @report_names, "deseq2_volcano_plot" );
      }
      for my $key ( keys %$pairs ) {
        push( @report_files, $deseq2taskname, $key . $suffix . "_DESeq2_sig.csv" );
        push( @report_names, "deseq2_" . $key );
      }
      push( @copy_files, $deseq2taskname, "_DESeq2.csv" );

      push( @copy_files, $deseq2taskname, "_DESeq2_sig.csv" );

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

      my $pairs = $config->{pairs};
      for my $key ( keys %$pairs ) {
        push( @report_files, $webgestaltTaskName, $key . "_geneontology_Biological_Process.txt" );
        push( @report_files, $webgestaltTaskName, $key . "_geneontology_Cellular_Component.txt" );
        push( @report_files, $webgestaltTaskName, $key . "_geneontology_Molecular_Function.txt" );
        push( @report_files, $webgestaltTaskName, $key . "_pathway_KEGG.txt" );
        push( @report_names, "WebGestalt_GO_BP_" . $key );
        push( @report_names, "WebGestalt_GO_CC_" . $key );
        push( @report_names, "WebGestalt_GO_MF_" . $key );
        push( @report_names, "WebGestalt_KEGG_" . $key );
      }
      $hasFunctionalEnrichment = 1;
    }

    if ( defined $gseaTaskName ) {
      push( @copy_files, $gseaTaskName, ".gsea\$" );

      my $pairs = $config->{pairs};
      for my $key ( keys %$pairs ) {
        push( @report_files, $gseaTaskName, $key . $suffix . ".*gsea.csv" );
        push( @report_names, "gsea_" . $key );
      }
      $hasFunctionalEnrichment = 1;
    }

    my $fcOptions= getValue( $def, "featureCount_option" );
    my $fcMultiMapping = ($fcOptions =~ /-m/) ? "TRUE" : "FALSE";
    my $options = {
      "DE_fold_change" => [ getValue( $def, "DE_fold_change", 2 ) ],
      "DE_pvalue"      => [ getValue( $def, "DE_pvalue",      0.05 ) ],
      "featureCounts_UseMultiMappingReads" => [$fcMultiMapping]
    };
    
    $config->{report} = {
      class                      => "CQS::BuildReport",
      perform                    => 1,
      target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . "report",
      report_rmd_file            => "../Pipeline/RNASeq.Rmd",
      additional_rmd_files       => "Functions.Rmd",
      parameterSampleFile1_ref   => \@report_files,
      parameterSampleFile1_names => \@report_names,
      parameterSampleFile2       => $options,
      parameterSampleFile3_ref   => \@copy_files,
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
