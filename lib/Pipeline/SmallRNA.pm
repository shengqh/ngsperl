#!/usr/bin/perl
package Pipeline::SmallRNA;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::SmallRNAUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performSmallRNA performSmallRNATask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.06';

sub getSmallRNAConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my ( $config, $individual_ref, $summary_ref, $cluster, $not_identical_ref, $preprocessing_dir, $class_independent_dir ) = getPrepareConfig( $def, 1 );

  my $host_genome_dir        = create_directory_or_die( $def->{target_dir} . "/host_genome" );
  my $nonhost_library_dir    = create_directory_or_die( $def->{target_dir} . "/nonhost_library" );
  my $nonhost_genome_dir     = create_directory_or_die( $def->{target_dir} . "/nonhost_genome" );
  my $nonhost_blast_dir      = create_directory_or_die( $def->{target_dir} . "/nonhost_blast" );
  my $data_visualization_dir = create_directory_or_die( $def->{target_dir} . "/data_visualization" );

  my @individual = @{$individual_ref};
  my @summary    = @{$summary_ref};

  my @table_for_correlation = ( "identical_sequence_count_table", ".count\$" );
  my @table_for_countSum    = ();
  my @table_for_readSummary    = ();
  my @table_for_pieSummary  = ( "identical", ".dupcount" );
  my @name_for_readSummary=();
  #print Dumper($config);

  my $search_not_identical  = ( !defined $def->{search_not_identical} )  || $def->{search_not_identical};
  my $search_host_genome    = ( !defined $def->{search_host_genome} )    || $def->{search_host_genome};
  my $search_miRBase        = ( !defined $def->{search_miRBase} )        || $def->{search_miRBase};
  my $search_unmapped_reads = ( !defined $def->{search_unmapped_reads} ) || $def->{search_unmapped_reads};
  my $blast_unmapped_reads = defined $def->{blast_unmapped_reads} && $def->{blast_unmapped_reads};
  my $do_comparison        = defined $def->{pairs};
  my $groups               = $def->{groups};
  my $groups_vis_layout    = $def->{groups_vis_layout};

  my $DE_show_gene_cluster        = $def->{DE_show_gene_cluster};
  my $DE_pvalue                   = $def->{DE_pvalue};
  my $DE_fold_change              = $def->{DE_fold_change};
  my $DE_add_count_one            = $def->{DE_add_count_one};
  my $DE_min_median_read_top100   = $def->{DE_min_median_read_top100};
  my $DE_min_median_read_smallRNA = $def->{DE_min_median_read_smallRNA};

  my $max_sequence_extension_base = $def->{max_sequence_extension_base};
  my $non_host_table_option       = "--maxExtensionBase $max_sequence_extension_base --outputReadTable --outputReadContigTable";

  my $blast_localdb = $def->{blast_localdb};

  if ($do_comparison) {
    my $class_independent = {
      deseq2_top100_reads => {
        class                => "Comparison::DESeq2",
        perform              => 1,
        target_dir           => $class_independent_dir . "/deseq2_top100_reads",
        option               => "",
        source_ref           => "pairs",
        groups_ref           => "groups",
        countfile_ref        => [ "identical_sequence_count_table", ".read.count\$" ],
        sh_direct            => 1,
        show_DE_gene_cluster => $DE_show_gene_cluster,
        pvalue               => $DE_pvalue,
        fold_change          => $DE_fold_change,
        min_median_read      => $DE_min_median_read_top100,
        add_count_one        => $DE_add_count_one,
        pbs                  => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "10",
          "mem"      => "10gb"
        },
      },
      deseq2_top100_contigs => {
        class                => "Comparison::DESeq2",
        perform              => 1,
        target_dir           => $class_independent_dir . "/deseq2_top100_contigs",
        option               => "",
        source_ref           => "pairs",
        groups_ref           => "groups",
        countfile_ref        => [ "identical_sequence_count_table", ".count\$" ],
        sh_direct            => 1,
        show_DE_gene_cluster => $DE_show_gene_cluster,
        pvalue               => $DE_pvalue,
        fold_change          => $DE_fold_change,
        min_median_read      => $DE_min_median_read_top100,
        add_count_one        => $DE_add_count_one,
        pbs                  => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "10",
          "mem"      => "10gb"
        },
      },
      deseq2_top100_reads_vis => {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $data_visualization_dir . "/deseq2_top100_reads_vis",
        rtemplate                => "DESeq2_all_vis.R",
        output_file              => ".Top100Reads.DESeq2.Matrix",
        output_file_ext          => ".png",
        parameterSampleFile1_ref => [ "deseq2_top100_reads", "_DESeq2.csv\$" ],
        parameterSampleFile2     => $def->{pairs_top100_deseq2_vis_layout},
        sh_direct                => 1,
        pbs                      => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },
      deseq2_top100_contigs_vis => {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $data_visualization_dir . "/deseq2_top100_contigs_vis",
        rtemplate                => "DESeq2_all_vis.R",
        output_file              => ".Top100Contigs.DESeq2.Matrix",
        output_file_ext          => ".png",
        parameterSampleFile1_ref => [ "deseq2_top100_contigs", "_DESeq2.csv\$" ],
        parameterSampleFile2     => $def->{pairs_top100_deseq2_vis_layout},
        sh_direct                => 1,
        pbs                      => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },
    };

    push @summary, ( "deseq2_top100_reads", "deseq2_top100_contigs", "deseq2_top100_reads_vis", "deseq2_top100_contigs_vis" );

    $config = merge( $config, $class_independent );
  }

  my $identical_ref = [ "identical", ".fastq.gz\$" ];

  if ($search_host_genome) {
    defined $def->{coordinate} or die "No smallRNA coordinate defined!";

    my $host_genome = {

      #1 mismatch search, NTA
      bowtie1_genome_1mm_NTA => {
        class         => "Alignment::Bowtie1",
        perform       => 1,
        target_dir    => $host_genome_dir . "/bowtie1_genome_1mm_NTA",
        option        => $def->{bowtie1_option_1mm},
        source_ref    => [ "identical_NTA", ".fastq.gz\$" ],
        bowtie1_index => $def->{bowtie1_index},
        samonly       => 0,
        sh_direct     => 1,
        mappedonly    => 1,
        cluster       => $cluster,
        pbs           => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=" . $def->{max_thread},
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      bowtie1_genome_1mm_NTA_smallRNA_count => {
        class           => "CQS::SmallRNACount",
        perform         => 1,
        target_dir      => $host_genome_dir . "/bowtie1_genome_1mm_NTA_smallRNA_count",
        option          => $def->{smallrnacount_option},
        source_ref      => "bowtie1_genome_1mm_NTA",
        fastq_files_ref => "identical_NTA",
        seqcount_ref    => [ "identical_NTA", ".dupcount\$" ],
        cqs_tools       => $def->{cqstools},
        coordinate_file => $def->{coordinate},
        fasta_file      => $def->{coordinate_fasta},
        sh_direct       => 1,
        cluster         => $cluster,
        pbs             => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      bowtie1_genome_1mm_NTA_smallRNA_table => {
        class      => "CQS::SmallRNATable",
        perform    => 1,
        target_dir => $host_genome_dir . "/bowtie1_genome_1mm_NTA_smallRNA_table",
        option     => "",
        source_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".mapped.xml" ],
        cqs_tools  => $def->{cqstools},
        prefix     => "smallRNA_1mm_",
        sh_direct  => 1,
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "10",
          "mem"      => "10gb"
        },
      },
      bowtie1_genome_1mm_NTA_smallRNA_category => {
        class                     => "CQS::UniqueR",
        perform                   => 1,
        target_dir                => $host_genome_dir . "/bowtie1_genome_1mm_NTA_smallRNA_category",
        rtemplate                 => "countTableVisFunctions.R,smallRnaCategory.R",
        output_file               => "",
        output_file_ext           => ".Category.Table.csv",
        parameterSampleFile1_ref  => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".info" ],
        parameterSampleFile2      => $groups,
        parameterSampleFile2Order => $def->{groups_order},
        parameterSampleFile3      => $groups_vis_layout,
        rCode                     => 'textSize=9;groupTextSize=' . $def->{table_vis_group_text_size} . ';',
        sh_direct                 => 1,
        pbs                       => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },
      host_genome_tRNA_category => {
        class                     => "CQS::UniqueR",
        perform                   => 1,
        target_dir                => $data_visualization_dir . "/host_genome_tRNA_category",
        rtemplate                 => "countTableVisFunctions.R,hostTrnaMappingVis.R",
        output_file               => ".tRNAMapping.Result",
        output_file_ext           => ".tRNAType2.Barplot.png",
        parameterSampleFile1Order => $def->{groups_order},
        parameterSampleFile1      => $groups,
        parameterSampleFile2      => $groups_vis_layout,
        parameterFile1_ref        => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.count\$" ],
        parameterFile3_ref        => [ "fastqc_count_vis", ".Reads.csv\$" ],
        rCode                     => 'maxCategory=3;textSize=9;groupTextSize=' . $def->{table_vis_group_text_size} . ';',
        sh_direct                 => 1,
        pbs                       => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },
    };

    push @table_for_pieSummary,  ( "bowtie1_genome_1mm_NTA_smallRNA_count", ".count\$" );
    push @table_for_correlation, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".count\$" );
    push @table_for_readSummary,
          ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.read.count\$", "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.read.count\$", "bowtie1_genome_1mm_NTA_smallRNA_table", ".other.read.count\$" );
    push @name_for_readSummary, ("Host miRNA","Host tRNA","Host other small RNA");
    push @table_for_countSum,
      ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.count\$", "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.count\$", "bowtie1_genome_1mm_NTA_smallRNA_table", ".other.count\$" );
    push @individual, ( "bowtie1_genome_1mm_NTA", "bowtie1_genome_1mm_NTA_smallRNA_count" );
    push @summary, ( "bowtie1_genome_1mm_NTA_smallRNA_table", "bowtie1_genome_1mm_NTA_smallRNA_category", "host_genome_tRNA_category" );

    if ($search_not_identical) {

      #not identical, for IGV
      $host_genome->{bowtie1_genome_1mm_notidentical} = {
        class         => "Alignment::Bowtie1",
        perform       => 1,
        target_dir    => $host_genome_dir . "/bowtie1_genome_1mm_notidentical",
        option        => $def->{bowtie1_option_1mm},
        source_ref    => $not_identical_ref,
        bowtie1_index => $def->{bowtie1_index},
        samonly       => 0,
        sh_direct     => 0,
        mappedonly    => 1,
        cluster       => $cluster,
        pbs           => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=" . $def->{max_thread},
          "walltime" => "72",
          "mem"      => "40gb"
        },
      };
      push @individual, ("bowtie1_genome_1mm_notidentical");
    }

    $config = merge( $config, $host_genome );

    if ($do_comparison) {
      my $comparison = {

        #DESeq2
        deseq2_miRNA => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $host_genome_dir . "/deseq2_miRNA",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
         deseq2_miRNA_NTA => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $host_genome_dir . "/deseq2_miRNA_NTA",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.NTA.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
         deseq2_miRNA_isomiR => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $host_genome_dir . "/deseq2_miRNA_isomiR",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.isomiR.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        deseq2_miRNA_isomiR_NTA => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $host_genome_dir . "/deseq2_miRNA_isomiR_NTA",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.isomiR_NTA.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        deseq2_tRNA => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $host_genome_dir . "/deseq2_tRNA",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        deseq2_tRNA_aminoacid => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $host_genome_dir . "/deseq2_tRNA_aminoacid",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.aminoacid.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        deseq2_otherSmallRNA => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $host_genome_dir . "/deseq2_otherSmallRNA",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".other.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        host_genome_deseq2_vis => {
          class                    => "CQS::UniqueR",
          perform                  => 1,
          target_dir               => $data_visualization_dir . "/host_genome_deseq2_vis",
          rtemplate                => "DESeq2_all_vis.R",
          output_file              => ".HostGenome.DESeq2.Matrix",
          output_file_ext          => ".png",
          parameterSampleFile1_ref => [ "deseq2_miRNA", "_DESeq2.csv\$", "deseq2_tRNA", "_DESeq2.csv\$", "deseq2_otherSmallRNA", "_DESeq2.csv\$" ],
          parameterSampleFile2     => $def->{pairs_host_deseq2_vis_layout},
          sh_direct                => 1,
          pbs                      => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "1",
            "mem"      => "10gb"
          },
        },
        host_genome_deseq2_miRNA_vis => {
          class                    => "CQS::UniqueR",
          perform                  => 1,
          target_dir               => $data_visualization_dir . "/host_genome_deseq2_miRNA_vis",
          rtemplate                => "DESeq2_all_vis.R",
          output_file              => ".HostGenome.miRNA.DESeq2.Matrix",
          output_file_ext          => ".png",
          parameterSampleFile1_ref => [ "deseq2_miRNA_isomiR", "_DESeq2.csv\$", "deseq2_miRNA_NTA", "_DESeq2.csv\$", "deseq2_miRNA_isomiR_NTA", "_DESeq2.csv\$" ],
          parameterSampleFile2     => $def->{pairs_host_deseq2_vis_layout},
          sh_direct                => 1,
          pbs                      => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "1",
            "mem"      => "10gb"
          },
        },
      };

      $config = merge( $config, $comparison );
      push @summary, ( "deseq2_miRNA", "deseq2_tRNA", "deseq2_tRNA_aminoacid", "deseq2_otherSmallRNA", "host_genome_deseq2_vis","deseq2_miRNA_NTA","deseq2_miRNA_isomiR","deseq2_miRNA_isomiR_NTA","host_genome_deseq2_miRNA_vis" );
    }

    if ( $do_comparison or defined $groups or defined $def->{tRNA_vis_group} ) {
      my $trna_vis_groups;
      my $trna_sig_result;
      if ( defined $def->{tRNA_vis_group} ) {
        $trna_vis_groups = $def->{tRNA_vis_group};
      }
      else {
        $trna_vis_groups = $groups;
      }
      if ($do_comparison) {
        $trna_sig_result = [ "deseq2_tRNA", "_DESeq2_sig.csv\$" ];
      }
      $config->{host_genome_tRNA_PositionVis} = {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $data_visualization_dir . "/host_genome_tRNA_PositionVis",
        rtemplate                => "tRNAPositionVis.R",
        output_file              => ".tRNAPositionVis",
        output_file_ext          => ".alltRNAPosition.png",
        parameterSampleFile1_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".tRNA.position\$" ],
        parameterSampleFile2     => $trna_vis_groups,
        parameterSampleFile3_ref => $trna_sig_result,
        sh_direct                => 1,
        pbs                      => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      };
      push @summary, ("host_genome_tRNA_PositionVis");
    }

    my $unmapped_reads = {

      #perfect matched reads with host genome
      bowtie1_genome_1mm_NTA_pmnames => {
        class      => "Samtools::PerfectMappedReadNames",
        perform    => 1,
        target_dir => $host_genome_dir . "/bowtie1_genome_1mm_NTA_pmnames",
        option     => "",
        source_ref => "bowtie1_genome_1mm_NTA",
        sh_direct  => 1,
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "10",
          "mem"      => "10gb"
        },
      },

      bowtie1_genome_unmapped_reads => {
        class            => "CQS::Perl",
        perform          => 1,
        target_dir       => $host_genome_dir . "/bowtie1_genome_unmapped_reads",
        perlFile         => "unmappedReadsToFastq.pl",
        source_ref       => [ "identical", ".fastq.gz\$" ],
        source2_ref      => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".mapped.xml" ],
        source3_ref      => ["bowtie1_genome_1mm_NTA_pmnames"],
        output_ext       => "_clipped_identical.unmapped.fastq.gz",
        output_other_ext => "_clipped_identical.unmapped.fastq.dupcount",
        sh_direct        => 1,
        pbs              => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      }
    };
    $config = merge( $config, $unmapped_reads );
    push @individual,           ( "bowtie1_genome_1mm_NTA_pmnames", "bowtie1_genome_unmapped_reads" );
    push @table_for_pieSummary, ( "bowtie1_genome_unmapped_reads",  ".dupcount" );
    $identical_ref = [ "bowtie1_genome_unmapped_reads", ".fastq.gz\$" ];
  }

  my @mapped  = ();
  my @pmnames = ();

  if ($search_miRBase) {
    my $mirbase = {
      bowtie1_miRBase_pm => {
        class         => "Alignment::Bowtie1",
        perform       => 1,
        target_dir    => $nonhost_library_dir . "/bowtie1_miRBase_pm",
        option        => $def->{bowtie1_option_pm},
        source_ref    => $identical_ref,
        bowtie1_index => $def->{bowtie1_miRBase_index},
        samonly       => 0,
        sh_direct     => 1,
        mappedonly    => 1,
        cluster       => $cluster,
        pbs           => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=" . $def->{max_thread},
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      bowtie1_miRBase_pm_count => {
        class        => "CQS::CQSChromosomeCount",
        perform      => 1,
        target_dir   => $nonhost_library_dir . "/bowtie1_miRBase_pm_count",
        option       => $def->{mirbase_count_option} . " -m --keepChrInName --keepSequence",
        source_ref   => "bowtie1_miRBase_pm",
        seqcount_ref => [ "identical", ".dupcount\$" ],
        cqs_tools    => $def->{cqstools},
        sh_direct    => 1,
        cluster      => $cluster,
        pbs          => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      bowtie1_miRBase_pm_table => {
        class      => "CQS::CQSChromosomeTable",
        perform    => 1,
        target_dir => $nonhost_library_dir . "/bowtie1_miRBase_pm_table",
        option     => $non_host_table_option,
        source_ref => [ "bowtie1_miRBase_pm_count", ".xml" ],
        cqs_tools  => $def->{cqstools},
        prefix     => "miRBase_pm_",
        sh_direct  => 1,
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "10",
          "mem"      => "10gb"
        },
      }
    };

    $config = merge( $config, $mirbase );

    #		push @table_for_correlation, ( "bowtie1_miRBase_pm_table", ".count\$" );
    push @table_for_countSum, ( "bowtie1_miRBase_pm_table", ".count\$" );
    push @individual,         ( "bowtie1_miRBase_pm",       "bowtie1_miRBase_pm_count" );
    push @summary,            ("bowtie1_miRBase_pm_table");

    push @mapped, ( "bowtie1_miRBase_pm_count", ".xml" );
  }

  if ($search_unmapped_reads) {
    my $unmappedreads = {

      # unmapped reads to tRNA
      bowtie1_tRNA_pm => {
        class         => 'Alignment::Bowtie1',
        cluster       => $cluster,
        sh_direct     => 1,
        perform       => 1,
        target_dir    => $nonhost_library_dir . "/bowtie1_tRNA_pm",
        samonly       => 0,
        source_ref    => $identical_ref,
        mappedonly    => 1,
        bowtie1_index => $def->{bowtie1_tRNA_index},
        option        => $def->{bowtie1_option_pm},
        pbs           => {
          'email'    => $def->{email},
          'walltime' => '72',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=8'
        }
      },

      bowtie1_tRNA_pm_count => {
        class        => 'CQS::CQSChromosomeCount',
        cluster      => $cluster,
        sh_direct    => 1,
        perform      => 1,
        target_dir   => $nonhost_library_dir . "/bowtie1_tRNA_pm_count",
        option       => $def->{smallrnacount_option} . " --keepChrInName --keepSequence",
        source_ref   => 'bowtie1_tRNA_pm',
        cqs_tools    => $def->{cqstools},
        seqcount_ref => [ "identical", ".dupcount\$" ],
        pbs          => {
          'email'    => $def->{email},
          'walltime' => '72',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=1'
        },
      },

      bowtie1_tRNA_pm_table => {
        class      => 'CQS::CQSChromosomeTable',
        cluster    => $cluster,
        sh_direct  => 1,
        perform    => 1,
        target_dir => $nonhost_library_dir . "/bowtie1_tRNA_pm_table",
        source_ref => [ 'bowtie1_tRNA_pm_count', '.xml' ],
        cqs_tools  => $def->{cqstools},
        option     => $non_host_table_option . ' --categoryMapFile ' . $def->{trna_category_map},
        prefix     => 'tRNA_pm_',
        pbs        => {
          'email'    => $def->{email},
          'walltime' => '10',
          'mem'      => '10gb',
          'nodes'    => '1:ppn=1'
        },
      },
      nonhost_library_tRNA_vis => {
        class                     => "CQS::UniqueR",
        perform                   => 1,
        target_dir                => $data_visualization_dir . "/nonhost_library_tRNA_vis",
        rtemplate                 => "countTableVisFunctions.R,bacteriaTrnaMappingVis.R",
        output_file               => ".tRNAMapping.Result",
        output_file_ext           => ".Species12.csv;.tRNAType1.csv;.tRNAType2.csv",
        parameterSampleFile1Order => $def->{groups_order},
        parameterSampleFile1      => $groups,
        parameterSampleFile2      => $groups_vis_layout,
        parameterFile1_ref        => [ "bowtie1_tRNA_pm_table", ".count\$" ],
        parameterFile3_ref        => [ "fastqc_count_vis", ".Reads.csv\$" ],
        rCode                     => 'maxCategory=3;textSize=9;groupTextSize=' . $def->{table_vis_group_text_size} . ';',
        sh_direct                 => 1,
        pbs                       => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },

      #unmapped reads to rRNA
      bowtie1_rRNA_pm => {
        pbs => {
          'email'    => $def->{email},
          'walltime' => '10',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=8'
        },
        cluster       => $cluster,
        sh_direct     => 1,
        perform       => 1,
        target_dir    => $nonhost_library_dir . "/bowtie1_rRNA_pm",
        samonly       => 0,
        mappedonly    => 1,
        source_ref    => $identical_ref,
        bowtie1_index => $def->{bowtie1_rRNA_index},
        option        => $def->{bowtie1_option_pm},
        class         => 'Alignment::Bowtie1'
      },
      bowtie1_rRNA_pm_count => {
        class        => 'CQS::CQSChromosomeCount',
        cluster      => $cluster,
        sh_direct    => 1,
        perform      => 1,
        target_dir   => $nonhost_library_dir . "/bowtie1_rRNA_pm_count",
        option       => $def->{smallrnacount_option} . ' --keepChrInName --keepSequence --categoryMapFile ' . $def->{rrna_category_map},
        source_ref   => 'bowtie1_rRNA_pm',
        cqs_tools    => $def->{cqstools},
        seqcount_ref => [ "identical", ".dupcount\$" ],
        pbs          => {
          'email'    => $def->{email},
          'walltime' => '72',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=1'
        },
      },

      bowtie1_rRNA_pm_table => {
        class      => 'CQS::CQSChromosomeTable',
        cluster    => $cluster,
        sh_direct  => 1,
        perform    => 1,
        target_dir => $nonhost_library_dir . "/bowtie1_rRNA_pm_table",
        source_ref => [ 'bowtie1_rRNA_pm_count', '.xml' ],
        cqs_tools  => $def->{cqstools},
        option     => $non_host_table_option,
        prefix     => 'rRNA_pm_',
        pbs        => {
          'email'    => $def->{email},
          'walltime' => '10',
          'mem'      => '10gb',
          'nodes'    => '1:ppn=1'
        },
      },
      nonhost_library_rRNA_vis => {
        class                     => "CQS::UniqueR",
        perform                   => 1,
        target_dir                => $data_visualization_dir . "/nonhost_library_rRNA_vis",
        rtemplate                 => "countTableVisFunctions.R,countTableVis.R",
        output_file               => ".rRNAMapping.Result",
        output_file_ext           => ".Barplot.png",
        parameterSampleFile1Order => $def->{groups_order},
        parameterSampleFile1      => $groups,
        parameterSampleFile2      => $groups_vis_layout,
        parameterFile1_ref        => [ "bowtie1_rRNA_pm_table", ".count\$" ],
        parameterFile3_ref        => [ "fastqc_count_vis", ".Reads.csv\$" ],
        rCode                     => 'maxCategory=NA;textSize=9;groupTextSize=' . $def->{table_vis_group_text_size} . ';',
        sh_direct                 => 1,
        pbs                       => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },

      #unmapped reads to group1 bacterial
      bowtie1_bacteria_group1_pm => {
        pbs => {
          'email'    => $def->{email},
          'walltime' => '72',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=8'
        },
        cluster       => $cluster,
        sh_direct     => 1,
        perform       => 1,
        target_dir    => $nonhost_genome_dir . "/bowtie1_bacteria_group1_pm",
        samonly       => 0,
        mappedonly    => 1,
        source_ref    => $identical_ref,
        bowtie1_index => $def->{bowtie1_bacteria_group1_index},
        option        => $def->{bowtie1_option_pm},
        class         => 'Alignment::Bowtie1'
      },

      bowtie1_bacteria_group1_pm_count => {
        pbs => {
          'email'    => $def->{email},
          'walltime' => '72',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=1'
        },
        cluster      => $cluster,
        sh_direct    => 1,
        perform      => 1,
        target_dir   => $nonhost_genome_dir . "/bowtie1_bacteria_group1_pm_count",
        option       => $def->{smallrnacount_option} . " --keepChrInName --keepSequence",
        source_ref   => 'bowtie1_bacteria_group1_pm',
        cqs_tools    => $def->{cqstools},
        seqcount_ref => [ "identical", ".dupcount\$" ],
        'class'      => 'CQS::CQSChromosomeCount'
      },

      bowtie1_bacteria_group1_pm_table => {
        class      => 'CQS::CQSChromosomeTable',
        cluster    => $cluster,
        sh_direct  => 1,
        perform    => 1,
        target_dir => $nonhost_genome_dir . "/bowtie1_bacteria_group1_pm_table",
        source_ref => [ 'bowtie1_bacteria_group1_pm_count', '.xml' ],
        cqs_tools  => $def->{cqstools},
        option     => $non_host_table_option . ' --categoryMapFile ' . $def->{bacteria_group1_species_map},
        prefix     => 'bacteria_group1_pm_',
        pbs        => {
          'email'    => $def->{email},
          'walltime' => '72',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=1'
        }
      },
      nonhost_genome_bacteria_group1_vis => {
        class                     => "CQS::UniqueR",
        perform                   => 1,
        target_dir                => $data_visualization_dir . "/nonhost_genome_bacteria_group1_vis",
        rtemplate                 => "countTableVisFunctions.R,countTableVis.R",
        output_file               => ".group1Mapping.Result",
        output_file_ext           => ".Piechart.png",
        parameterSampleFile1Order => $def->{groups_order},
        parameterSampleFile1      => $groups,
        parameterSampleFile2      => $groups_vis_layout,
        parameterFile1_ref        => [ "bowtie1_bacteria_group1_pm_table", ".category.count\$" ],
#        parameterFile2            => $def->{bacteria_group1_log},
        parameterFile3_ref        => [ "fastqc_count_vis", ".Reads.csv\$" ],
        rCode                     => 'maxCategory=4;textSize=9;groupTextSize=' . $def->{table_vis_group_text_size} . ';',
        sh_direct                 => 1,
        pbs                       => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },

      #unmapped reads to group2 bacterial
      bowtie1_bacteria_group2_pm => {
        pbs => {
          'email'    => $def->{email},
          'walltime' => '72',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=8'
        },
        cluster       => $cluster,
        sh_direct     => 1,
        perform       => 1,
        target_dir    => $nonhost_genome_dir . "/bowtie1_bacteria_group2_pm",
        samonly       => 0,
        mappedonly    => 1,
        source_ref    => $identical_ref,
        bowtie1_index => $def->{bowtie1_bacteria_group2_index},
        option        => $def->{bowtie1_option_pm},
        class         => 'Alignment::Bowtie1'
      },

      bowtie1_bacteria_group2_pm_count => {
        pbs => {
          'email'    => $def->{email},
          'walltime' => '72',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=1'
        },
        cluster      => $cluster,
        sh_direct    => 1,
        perform      => 1,
        target_dir   => $nonhost_genome_dir . "/bowtie1_bacteria_group2_pm_count",
        option       => $def->{smallrnacount_option} . " --keepChrInName --keepSequence",
        source_ref   => 'bowtie1_bacteria_group2_pm',
        cqs_tools    => $def->{cqstools},
        seqcount_ref => [ "identical", ".dupcount\$" ],
        'class'      => 'CQS::CQSChromosomeCount'
      },

      bowtie1_bacteria_group2_pm_table => {
        class      => 'CQS::CQSChromosomeTable',
        cluster    => $cluster,
        sh_direct  => 1,
        perform    => 1,
        target_dir => $nonhost_genome_dir . "/bowtie1_bacteria_group2_pm_table",
        source_ref => [ 'bowtie1_bacteria_group2_pm_count', '.xml' ],
        cqs_tools  => $def->{cqstools},
        option     => $non_host_table_option . ' --categoryMapFile ' . $def->{bacteria_group2_species_map},
        prefix     => 'bacteria_group2_pm_',
        pbs        => {
          'email'    => $def->{email},
          'walltime' => '72',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=1'
        }
      },

      nonhost_genome_bacteria_group2_vis => {
        class                     => "CQS::UniqueR",
        perform                   => 1,
        target_dir                => $data_visualization_dir . "/nonhost_genome_bacteria_group2_vis",
        rtemplate                 => "countTableVisFunctions.R,countTableVis.R",
        output_file               => ".group2Mapping.Result",
        output_file_ext           => ".Piechart.png",
        parameterSampleFile1Order => $def->{groups_order},
        parameterSampleFile1      => $groups,
        parameterSampleFile2      => $groups_vis_layout,
        parameterFile1_ref        => [ "bowtie1_bacteria_group2_pm_table", ".category.count\$" ],
#        parameterFile2            => $def->{bacteria_group2_log},
        parameterFile3_ref        => [ "fastqc_count_vis", ".Reads.csv\$" ],
        rCode                     => 'maxCategory=5;textSize=9;groupTextSize=' . $def->{table_vis_group_text_size} . ';',
        sh_direct                 => 1,
        pbs                       => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },

      #unmapped reads to group4 fungus
      bowtie1_fungus_group4_pm => {
        pbs => {
          'email'    => $def->{email},
          'walltime' => '72',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=8'
        },
        cluster       => $cluster,
        sh_direct     => 1,
        perform       => 1,
        target_dir    => $nonhost_genome_dir . "/bowtie1_fungus_group4_pm",
        samonly       => 0,
        mappedonly    => 1,
        source_ref    => $identical_ref,
        bowtie1_index => $def->{bowtie1_fungus_group4_index},
        option        => $def->{bowtie1_option_pm},
        class         => 'Alignment::Bowtie1'
      },

      bowtie1_fungus_group4_pm_count => {
        pbs => {
          'email'    => $def->{email},
          'walltime' => '72',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=1'
        },
        cluster      => $cluster,
        sh_direct    => 1,
        perform      => 1,
        target_dir   => $nonhost_genome_dir . "/bowtie1_fungus_group4_pm_count",
        option       => $def->{smallrnacount_option} . " --keepChrInName --keepSequence",
        source_ref   => 'bowtie1_fungus_group4_pm',
        cqs_tools    => $def->{cqstools},
        seqcount_ref => [ "identical", ".dupcount\$" ],
        'class'      => 'CQS::CQSChromosomeCount'
      },

      bowtie1_fungus_group4_pm_table => {
        class      => 'CQS::CQSChromosomeTable',
        cluster    => $cluster,
        sh_direct  => 1,
        perform    => 1,
        target_dir => $nonhost_genome_dir . "/bowtie1_fungus_group4_pm_table",
        source_ref => [ 'bowtie1_fungus_group4_pm_count', '.xml' ],
        cqs_tools  => $def->{cqstools},
        option     => $non_host_table_option . ' --categoryMapFile ' . $def->{fungus_group4_species_map},
        prefix     => 'fungus_group4_pm_',
        pbs        => {
          'email'    => $def->{email},
          'walltime' => '72',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=1'
        }
      },
      nonhost_genome_fungus_group4_vis => {
        class                     => "CQS::UniqueR",
        perform                   => 1,
        target_dir                => $data_visualization_dir . "/nonhost_genome_fungus_group4_vis",
        rtemplate                 => "countTableVisFunctions.R,countTableVis.R",
        output_file               => ".group4Mapping.Result",
        output_file_ext           => ".Piechart.png",
        parameterSampleFile1Order => $def->{groups_order},
        parameterSampleFile1      => $groups,
        parameterSampleFile2      => $groups_vis_layout,
        parameterFile1_ref        => [ "bowtie1_fungus_group4_pm_table", ".category.count\$" ],
#        parameterFile2            => $def->{fungus_group4_log},
        parameterFile3_ref        => [ "fastqc_count_vis", ".Reads.csv\$" ],
        sh_direct                 => 1,
        rCode                     => 'maxCategory=8;textSize=9;groupTextSize=' . $def->{table_vis_group_text_size} . ';',
        pbs                       => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },
      nonhost_overlap_vis => {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $data_visualization_dir . "/nonhost_overlap_vis",
        rtemplate                => "countTableVisFunctions.R,NonHostOverlap.R",
        output_file              => ".NonHost.Reads",
        output_file_ext          => ".Overlap.csv",
        parameterSampleFile1_ref => [
          "bowtie1_bacteria_group1_pm_table", ".read.count\$", "bowtie1_bacteria_group2_pm_table", ".read.count\$", "bowtie1_fungus_group4_pm_table", ".read.count\$",
          "bowtie1_tRNA_pm_table",            ".read.count\$", "bowtie1_rRNA_pm_table",            ".read.count\$",
        ],
        parameterSampleFile2Order => $def->{groups_order},
        parameterSampleFile2      => $groups,
        parameterSampleFile3      => $groups_vis_layout,

        #				parameterFile1_ref        => [ "bowtie1_fungus_group4_pm_table", ".count\$" ],
        #				parameterFile2            => $def->{fungus_group4_log},
        parameterFile3_ref => [ "fastqc_count_vis", ".Reads.csv\$" ],
        sh_direct          => 1,
        rCode              => 'maxCategory=8;textSize=9;groupTextSize=' . $def->{table_vis_group_text_size} . ';',
        pbs                => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },
    };

    $config = merge( $config, $unmappedreads );

    push @table_for_correlation,
      (
      "bowtie1_tRNA_pm_table",              ".count\$",       "bowtie1_rRNA_pm_table",            ".count\$", "nonhost_genome_bacteria_group1_vis", ".Species.csv\$",
      "nonhost_genome_bacteria_group2_vis", ".Species.csv\$", "nonhost_genome_fungus_group4_vis", ".Species.csv\$"
      );
    push @table_for_countSum,
      (
      "bowtie1_tRNA_pm_table",            ".category.count\$", "bowtie1_rRNA_pm_table",          ".count\$", "bowtie1_bacteria_group1_pm_table", ".count\$",
      "bowtie1_bacteria_group2_pm_table", ".count\$",          "bowtie1_fungus_group4_pm_table", ".count\$"
      );
    push @table_for_readSummary,      
    (
      "bowtie1_tRNA_pm_table",              ".read.count\$",       "bowtie1_rRNA_pm_table",            ".read.count\$", "bowtie1_bacteria_group1_pm_table", ".read.count\$",
      "bowtie1_bacteria_group2_pm_table", ".read.count\$", "bowtie1_fungus_group4_pm_table", ".read.count\$"
      );
    push @name_for_readSummary, ("Non host tRNA","Non host rRNA","Human Microbiome Bacteria","Environment Bacteria","Fungus");
 
    push @individual,
      (
      "bowtie1_tRNA_pm",            "bowtie1_tRNA_pm_count",            "bowtie1_rRNA_pm",            "bowtie1_rRNA_pm_count",
      "bowtie1_bacteria_group1_pm", "bowtie1_bacteria_group1_pm_count", "bowtie1_bacteria_group2_pm", "bowtie1_bacteria_group2_pm_count",
      "bowtie1_fungus_group4_pm",   "bowtie1_fungus_group4_pm_count"
      );
    push @summary,
      (
      "bowtie1_tRNA_pm_table",            "nonhost_library_tRNA_vis",           "bowtie1_rRNA_pm_table",            "nonhost_library_rRNA_vis",
      "bowtie1_bacteria_group1_pm_table", "nonhost_genome_bacteria_group1_vis", "bowtie1_bacteria_group2_pm_table", "nonhost_genome_bacteria_group2_vis",
      "bowtie1_fungus_group4_pm_table",   "nonhost_genome_fungus_group4_vis",   "nonhost_overlap_vis"
      );

    push @mapped,
      (
      "bowtie1_tRNA_pm_count",            ".xml", "bowtie1_rRNA_pm_count",          ".xml", "bowtie1_bacteria_group1_pm_count", ".xml",
      "bowtie1_bacteria_group2_pm_count", ".xml", "bowtie1_fungus_group4_pm_count", ".xml"
      );

    #do unmapped reads DESeq2
    if ($do_comparison) {
      my $unmapped_comparison = {

        #DESeq2
        deseq2_nonhost_tRNA => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $nonhost_library_dir . "/deseq2_nonhost_tRNA",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_tRNA_pm_table", ".count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        deseq2_nonhost_tRNA_category => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $nonhost_library_dir . "/deseq2_nonhost_tRNA_category",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_tRNA_pm_table", ".category.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        deseq2_nonhost_tRNA_species => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $nonhost_library_dir . "/deseq2_nonhost_tRNA_species",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "nonhost_library_tRNA_vis", ".Species12.csv\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        deseq2_nonhost_tRNA_type => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $nonhost_library_dir . "/deseq2_nonhost_tRNA_type",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "nonhost_library_tRNA_vis", ".tRNAType1.csv\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        deseq2_nonhost_tRNA_anticodon => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $nonhost_library_dir . "/deseq2_nonhost_tRNA_anticodon",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "nonhost_library_tRNA_vis", ".tRNAType2.csv\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        deseq2_nonhost_rRNA => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $nonhost_library_dir . "/deseq2_nonhost_rRNA",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_rRNA_pm_table", ".count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        nonhost_library_deseq2_vis => {
          class      => "CQS::UniqueR",
          perform    => 1,
          target_dir => $data_visualization_dir . "/nonhost_library_deseq2_vis",
          ,
          rtemplate                => "DESeq2_all_vis.R",
          output_file              => ".NonHostLibrary.DESeq2.Matrix",
          output_file_ext          => ".png",
          parameterSampleFile1_ref => [ "deseq2_nonhost_tRNA", "_DESeq2.csv\$", "deseq2_nonhost_rRNA", "_DESeq2.csv\$" ],
          parameterSampleFile2     => $def->{pairs_nonHosttRNArRNA_deseq2_vis_layout},
          sh_direct                => 1,
          pbs                      => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "1",
            "mem"      => "10gb"
          },
        },

        deseq2_bacteria_group1 => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $nonhost_genome_dir . "/deseq2_bacteria_group1",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_bacteria_group1_pm_table", ".category.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        deseq2_bacteria_group1_reads => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $nonhost_genome_dir . "/deseq2_bacteria_group1_reads",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_bacteria_group1_pm_table", ".read.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        deseq2_bacteria_group2 => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $nonhost_genome_dir . "/deseq2_bacteria_group2",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_bacteria_group2_pm_table", ".category.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        deseq2_bacteria_group2_reads => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $nonhost_genome_dir . "/deseq2_bacteria_group2_reads",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_bacteria_group2_pm_table", ".read.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        deseq2_fungus_group4 => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $nonhost_genome_dir . "/deseq2_fungus_group4",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_fungus_group4_pm_table", ".category.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        deseq2_fungus_group4_reads => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $nonhost_genome_dir . "/deseq2_fungus_group4_reads",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_fungus_group4_pm_table", ".read.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => $DE_show_gene_cluster,
          pvalue               => $DE_pvalue,
          fold_change          => $DE_fold_change,
          min_median_read      => $DE_min_median_read_smallRNA,
          add_count_one        => $DE_add_count_one,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        nonhost_genome_deseq2_vis => {
          class                    => "CQS::UniqueR",
          perform                  => 1,
          target_dir               => $data_visualization_dir . "/nonhost_genome_deseq2_vis",
          rtemplate                => "DESeq2_all_vis.R",
          output_file              => ".NonHostGenome.DESeq2.Matrix",
          output_file_ext          => ".png",
          parameterSampleFile1_ref => [ "deseq2_bacteria_group1", "_DESeq2.csv\$", "deseq2_bacteria_group2", "_DESeq2.csv\$", "deseq2_fungus_group4", "_DESeq2.csv\$" ],
          parameterSampleFile2     => $def->{pairs_nonHostGroups_deseq2_vis_layout},
          sh_direct                => 1,
          pbs                      => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "1",
            "mem"      => "10gb"
          },
        },
        nonhost_genome_deseq2_reads_vis => {
          class                    => "CQS::UniqueR",
          perform                  => 1,
          target_dir               => $data_visualization_dir . "/nonhost_genome_deseq2_reads_vis",
          rtemplate                => "DESeq2_all_vis.R",
          output_file              => ".NonHostGenomeReads.DESeq2.Matrix",
          output_file_ext          => ".png",
          parameterSampleFile1_ref => [ "deseq2_bacteria_group1_reads", "_DESeq2.csv\$", "deseq2_bacteria_group2_reads", "_DESeq2.csv\$", "deseq2_fungus_group4_reads", "_DESeq2.csv\$" ],
          parameterSampleFile2     => $def->{pairs_nonHostGroups_deseq2_vis_layout},
          sh_direct                => 1,
          pbs                      => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "1",
            "mem"      => "10gb"
          },
        },
      };

      $config = merge( $config, $unmapped_comparison );
      push @summary,
        (
        "deseq2_nonhost_tRNA",        "deseq2_nonhost_tRNA_species", "deseq2_nonhost_tRNA_type",     "deseq2_nonhost_tRNA_anticodon",
        "deseq2_nonhost_rRNA",        "nonhost_library_deseq2_vis",  "deseq2_bacteria_group1",       "deseq2_bacteria_group2",
        "deseq2_fungus_group4",       "nonhost_genome_deseq2_vis",   "deseq2_bacteria_group1_reads", "deseq2_bacteria_group2_reads",
        "deseq2_fungus_group4_reads", "nonhost_genome_deseq2_reads_vis"
        );
    }
  }

  if ( $search_miRBase || $search_unmapped_reads ) {
    my $unmapped_reads = {
      bowtie1_unmapped_reads => {
        class            => "CQS::Perl",
        perform          => 1,
        target_dir       => $nonhost_blast_dir . "/bowtie1_unmapped_reads",
        perlFile         => "unmappedReadsToFastq.pl",
        source_ref       => $identical_ref,
        source2_ref      => \@mapped,
        source3_ref      => \@pmnames,
        output_ext       => "_clipped_identical.unmapped.fastq.gz",
        output_other_ext => "_clipped_identical.unmapped.fastq.dupcount",
        sh_direct        => 1,
        pbs              => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },
    };

    $identical_ref = [ "bowtie1_unmapped_reads", ".fastq.gz\$" ];
    $config = merge( $config, $unmapped_reads );
    push @individual, ("bowtie1_unmapped_reads");
    push @table_for_pieSummary, ( "bowtie1_unmapped_reads", ".dupcount" );
  }

  if ($blast_unmapped_reads) {
    my $blast = {

      bowtie1_unmapped_sequence_count_table => {
        class           => "CQS::SmallRNASequenceCountTable",
        perform         => 1,
        target_dir      => $nonhost_blast_dir . "/bowtie1_unmapped_sequence_count_table",
        option          => "",
        source_ref      => [ "identical", ".dupcount\$" ],
        fastq_files_ref => $identical_ref,
        cqs_tools       => $def->{cqstools},
        suffix          => "_unmapped",
        sh_direct       => 1,
        cluster         => $cluster,
        pbs             => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "10",
          "mem"      => "10gb"
        },
      },
      bowtie1_unmapped_sequence_blast => {
        class      => "Blast::Blastn",
        perform    => 1,
        target_dir => $nonhost_blast_dir . "/bowtie1_unmapped_sequence_blast",
        option     => "",
        source_ref => [ "bowtie1_unmapped_sequence_count_table", ".fasta\$" ],
        sh_direct  => 0,
        localdb    => $blast_localdb,
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=" . $def->{max_thread},
          "walltime" => "10",
          "mem"      => "10gb"
        },
      },
    };

    $config = merge( $config, $blast );
    push @summary, ( "bowtie1_unmapped_sequence_count_table", "bowtie1_unmapped_sequence_blast" );
  }

  $config->{count_table_correlation} = {
    class                     => "CQS::UniqueR",
    perform                   => 1,
    target_dir                => $data_visualization_dir . "/count_table_correlation",
    rtemplate                 => "countTableVisFunctions.R,countTableGroupCorrelation.R",
    output_file               => "parameterSampleFile1",
    output_file_ext           => ".Correlation.png",
    parameterSampleFile1_ref  => \@table_for_correlation,
    parameterSampleFile2Order => $def->{groups_order},
    parameterSampleFile2      => $groups,
    sh_direct                 => 1,
    pbs                       => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "10gb"
    },
    },
    $config->{reads_in_tasks} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $data_visualization_dir . "/reads_in_tasks",
    rtemplate                => "countTableVisFunctions.R,ReadsInTasks.R",
    output_file_ext          => ".TaskReads.csv",
    parameterSampleFile1_ref => \@table_for_countSum,
    parameterFile3_ref       => [ "fastqc_count_vis", ".Reads.csv\$" ],
    sh_direct                => 1,
    pbs                      => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "12",
      "mem"      => "10gb"
    },
    };
  $config->{reads_in_tasks_pie} = {
    class                    => "CQS::UniqueR",
    suffix                   => "_pie",
    perform                  => 1,
    target_dir               => $data_visualization_dir . "/reads_in_tasks",
    rtemplate                => "countTableVisFunctions.R,ReadsInTasksPie.R",
    output_file_ext          => ".NonParallel.TaskReads.csv",
    parameterSampleFile1_ref => \@table_for_pieSummary,
    parameterSampleFile2      => $groups,
    parameterSampleFile3      => $groups_vis_layout,
    #    parameterFile3_ref       => [ "fastqc_count_vis", ".Reads.csv\$" ],
    sh_direct => 1,
    pbs       => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "12",
      "mem"      => "10gb"
    },
  };
  my $name_for_readSummary_r="readFilesModule=c('".join("','",@name_for_readSummary)."')";
  $config->{reads_mapping_summary} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $data_visualization_dir . "/reads_mapping_summary",
    rtemplate                => "countTableVisFunctions.R,ReadsMappingSummary.R",
    output_file_ext          => ".ReadsMapping.Summary.csv",
    parameterFile1_ref       =>[ "identical_sequence_count_table", ".read.count\$" ],
    parameterSampleFile1_ref => \@table_for_readSummary,
    parameterSampleFile2      => $groups,
    parameterSampleFile3      => $groups_vis_layout,
    rCode=>$name_for_readSummary_r,
    sh_direct => 1,
    pbs       => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "12",
      "mem"      => "10gb"
    },
  };
  push @summary, ( "count_table_correlation", "reads_in_tasks", "reads_in_tasks_pie","reads_mapping_summary" );

  $config->{sequencetask} = {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => $def->{target_dir} . "/sequencetask",
    option     => "",
    source     => {
      step1 => \@individual,
      step2 => \@summary,
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

sub performSmallRNA {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getSmallRNAConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

sub performSmallRNATask {
  my ( $def, $task ) = @_;

  my $config = getSmallRNAConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
