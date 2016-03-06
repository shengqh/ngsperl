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

our $VERSION = '0.02';

sub getSmallRNAConfig {
  my ($def) = @_;
  $def->{VERSION}=$VERSION;
  
  my ( $config, $individual_ref, $summary_ref, $cluster, $not_identical_ref ) = getPrepareConfig( $def, 1 );
  my @individual = @{$individual_ref};
  my @summary    = @{$summary_ref};

  #print Dumper($config);

  my $search_not_identical  = ( !defined $def->{search_not_identical} ) || $def->{search_not_identical};
  my $search_host_genome    = defined $def->{bowtie1_index};
  my $search_miRBase        = defined $def->{bowtie1_miRBase_index};
  my $search_unmapped_reads = ( !defined $def->{search_unmapped_reads} ) || $def->{search_unmapped_reads};
  my $blast_unmapped_reads  = defined $def->{blast_unmapped_reads} && $def->{blast_unmapped_reads};
  my $do_comparison         = defined $def->{pairs};
  my $groups                = $def->{groups};

  if ($do_comparison) {
    $config->{top100Reads_deseq2} = {
      class                => "Comparison::DESeq2",
      perform              => 1,
      target_dir           => $def->{target_dir} . "/top100Reads_deseq2",
      option               => "",
      source_ref           => "pairs",
      groups_ref           => "groups",
      countfile_ref        => [ "identical_sequence_count_table", ".count\$" ],
      sh_direct            => 1,
      show_DE_gene_cluster => 1,
      pvalue               => 0.05,
      fold_change          => 1.5,
      min_median_read      => 1,
      pbs                  => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    };
    push @summary, ("top100Reads_deseq2");
  }

  my $identical_ref = [ "identical", ".fastq.gz\$" ];

  if ($search_host_genome) {
    defined $def->{coordinate} or die "No smallRNA coordinate defined!";

    my $host_genome = {

      #1 mismatch search, NTA
      bowtie1_genome_1mm_NTA => {
        class         => "Alignment::Bowtie1",
        perform       => 1,
        target_dir    => $def->{target_dir} . "/bowtie1_genome_1mm_NTA",
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
        target_dir      => $def->{target_dir} . "/bowtie1_genome_1mm_NTA_smallRNA_count",
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
        target_dir => $def->{target_dir} . "/bowtie1_genome_1mm_NTA_smallRNA_table",
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
      bowtie1_genome_1mm_NTA_smallRNA_table_vis => {
        class                => "CQS::UniqueR",
        perform              => 1,
        target_dir           => $def->{target_dir} . "/bowtie1_genome_1mm_NTA_smallRNA_table",
        rtemplate            => "countTableCorrelation.R",
        output_file          => "parameterSampleFile1",
        output_file_ext      => ".Correlation.png",
        parameterSampleFile1_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".count\$" ],
        sh_direct            => 1,
        pbs                  => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },
      bowtie1_genome_1mm_NTA_smallRNA_category => {
        class      => "CQS::SmallRNACategory",
        perform    => 1,
        target_dir => $def->{target_dir} . "/bowtie1_genome_1mm_NTA_smallRNA_category",
        option     => "",
        source_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".info\$" ],
        cqs_tools  => $def->{cqstools},
        sh_direct  => 1,
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
    };
    push @individual, ( "bowtie1_genome_1mm_NTA", "bowtie1_genome_1mm_NTA_smallRNA_count" );
    push @summary, ( "bowtie1_genome_1mm_NTA_smallRNA_table", "bowtie1_genome_1mm_NTA_smallRNA_table_vis","bowtie1_genome_1mm_NTA_smallRNA_category" );

    if ($search_not_identical) {

      #not identical, for IGV
      $host_genome->{bowtie1_genome_1mm_notidentical} = {
        class         => "Alignment::Bowtie1",
        perform       => 1,
        target_dir    => $def->{target_dir} . "/bowtie1_genome_1mm_notidentical",
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
        miRNA_deseq2 => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $def->{target_dir} . "/miRNA_deseq2",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => 1,
          pvalue               => 0.05,
          fold_change          => 1.5,
          min_median_read      => 5,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        tRNA_deseq2 => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $def->{target_dir} . "/tRNA_deseq2",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => 1,
          pvalue               => 0.05,
          fold_change          => 1.5,
          min_median_read      => 5,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        otherSmallRNA_deseq2 => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $def->{target_dir} . "/otherSmallRNA_deseq2",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".other.count\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => 1,
          pvalue               => 0.05,
          fold_change          => 1.5,
          min_median_read      => 5,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
      };

      $config = merge( $config, $comparison );
      push @summary, ( "miRNA_deseq2", "tRNA_deseq2", "otherSmallRNA_deseq2" );
    }

    if ( $do_comparison or defined $def->{tRNA_vis_group} ) {
      my $trna_vis_groups;
      my $trna_sig_result;
      if ( defined $def->{tRNA_vis_group} ) {
        $trna_vis_groups = $def->{tRNA_vis_group};
      }
      else {
        $trna_vis_groups = $def->{groups};
      }
      if ($do_comparison) {
        $trna_sig_result = [ "tRNA_deseq2", "_DESeq2_sig.csv\$" ];
      }
      $config->{tRNA_PositionVis} = {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $def->{target_dir} . "/tRNA_PositionVis",
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
      push @summary, ("tRNA_PositionVis");
    }

    # extract unmapped reads. the reads mapped to host smallRNA with/without mismatch and perfect mapped to host genome will be excluded.
    if ( $search_miRBase || $search_unmapped_reads || $blast_unmapped_reads ) {
      my $unmapped_reads = {

        #perfect matched reads with host genome
        bowtie1_genome_1mm_NTA_pmnames => {
          class      => "Samtools::PerfectMappedReadNames",
          perform    => 1,
          target_dir => $def->{target_dir} . "/bowtie1_genome_1mm_NTA_pmnames",
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
          class       => "CQS::Perl",
          perform     => 1,
          target_dir  => $def->{target_dir} . "/bowtie1_genome_unmapped_reads",
          perlFile    => "unmappedReadsToFastq.pl",
          source_ref  => [ "identical", ".fastq.gz\$" ],
          source2_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".mapped.xml" ],
          source3_ref => ["bowtie1_genome_1mm_NTA_pmnames"],
          output_ext  => "_clipped_identical.unmapped.fastq.gz",
          output_other_ext  => "_clipped_identical.unmapped.fastq.dupcount",
          sh_direct   => 1,
          pbs         => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "1",
            "mem"      => "10gb"
          },
        }
      };
      $config = merge( $config, $unmapped_reads );
      push @individual, ( "bowtie1_genome_1mm_NTA_pmnames", "bowtie1_genome_unmapped_reads" );
      $identical_ref = [ "bowtie1_genome_unmapped_reads", ".fastq.gz\$" ];
    }
  }

  my @mapped  = ();
  my @pmnames = ();

  if ($search_miRBase) {
    my $mirbase = {
      bowtie1_miRBase_pm => {
        class         => "Alignment::Bowtie1",
        perform       => 1,
        target_dir    => $def->{target_dir} . "/bowtie1_miRBase_pm",
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
        target_dir   => $def->{target_dir} . "/bowtie1_miRBase_pm_count",
        option       => $def->{mirbase_count_option} . " -m",
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
        target_dir => $def->{target_dir} . "/bowtie1_miRBase_pm_table",
        option     => "",
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
    push @individual, ( "bowtie1_miRBase_pm", "bowtie1_miRBase_pm_count" );
    push @summary, ("bowtie1_miRBase_pm_table");

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
        target_dir    => $def->{target_dir} . "/bowtie1_tRNA_pm",
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
        target_dir   => $def->{target_dir} . "/bowtie1_tRNA_pm_count",
        option       => $def->{smallrnacount_option},
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
        target_dir => $def->{target_dir} . "/bowtie1_tRNA_pm_table",
        source_ref => [ 'bowtie1_tRNA_pm_count', '.xml' ],
        cqs_tools  => $def->{cqstools},
        option     => '',
        prefix     => 'tRNA_pm_',
        pbs        => {
          'email'    => $def->{email},
          'walltime' => '10',
          'mem'      => '10gb',
          'nodes'    => '1:ppn=1'
        },
      },
      bowtie1_tRNA_pm_table_vis => {
        class                => "CQS::UniqueR",
        perform              => 1,
        target_dir           => $def->{target_dir} . "/bowtie1_tRNA_pm_table",
        rtemplate            => "BacTrnaMappingVis.R",
        output_file          => "",
        output_file_ext      => ".top.png",
        parameterSampleFile1 => $groups,
        parameterFile1_ref   => [ "bowtie1_tRNA_pm_table", ".count\$" ],
        sh_direct            => 1,
        pbs                  => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },

      #unmapped reads to rRNAL
      bowtie1_rRNAL_pm => {
        pbs => {
          'email'    => $def->{email},
          'walltime' => '10',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=8'
        },
        cluster       => $cluster,
        sh_direct     => 1,
        perform       => 1,
        target_dir    => $def->{target_dir} . "/bowtie1_rRNAL_pm",
        samonly       => 0,
        mappedonly    => 1,
        source_ref    => $identical_ref,
        bowtie1_index => $def->{bowtie1_rRNAL_index},
        option        => $def->{bowtie1_option_pm},
        class         => 'Alignment::Bowtie1'
      },

      bowtie1_rRNAL_pm_names => {
        class      => "Samtools::PerfectMappedReadNames",
        perform    => 1,
        target_dir => $def->{target_dir} . "/bowtie1_rRNAL_pm_names",
        option     => "",
        source_ref => "bowtie1_rRNAL_pm",
        sh_direct  => 1,
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "10",
          "mem"      => "10gb"
        },
      },

      bowtie1_rRNAL_pm_table => {
        pbs => {
          'email'    => $def->{email},
          'walltime' => '10',
          'mem'      => '10gb',
          'nodes'    => '1:ppn=1'
        },
        cluster      => $cluster,
        sh_direct    => 1,
        perform      => 1,
        target_dir   => $def->{target_dir} . "/bowtie1_rRNAL_pm_table",
        source_ref   => 'bowtie1_rRNAL_pm',
        seqcount_ref => [ "identical", ".dupcount\$" ],
        cqs_tools    => $def->{cqstools},
        option       => '',
        class        => 'CQS::BAMSequenceCountTable',
        prefix       => 'rRNAL_pm_'
      },

      #unmapped reads to rRNAS
      bowtie1_rRNAS_pm => {
        pbs => {
          'email'    => $def->{email},
          'walltime' => '10',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=8'
        },
        cluster       => $cluster,
        sh_direct     => 1,
        perform       => 1,
        target_dir    => $def->{target_dir} . "/bowtie1_rRNAS_pm",
        samonly       => 0,
        mappedonly    => 1,
        source_ref    => $identical_ref,
        bowtie1_index => $def->{bowtie1_rRNAS_index},
        option        => $def->{bowtie1_option_pm},
        class         => 'Alignment::Bowtie1'
      },
      bowtie1_rRNAS_pm_names => {
        class      => "Samtools::PerfectMappedReadNames",
        perform    => 1,
        target_dir => $def->{target_dir} . "/bowtie1_rRNAS_pm_names",
        option     => "",
        source_ref => "bowtie1_rRNAS_pm",
        sh_direct  => 1,
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "10",
          "mem"      => "10gb"
        },
      },

      bowtie1_rRNAS_pm_table => {
        pbs => {
          'email'    => $def->{email},
          'walltime' => '10',
          'mem'      => '10gb',
          'nodes'    => '1:ppn=1'
        },
        cluster      => $cluster,
        sh_direct    => 1,
        perform      => 1,
        target_dir   => $def->{target_dir} . "/bowtie1_rRNAS_pm_table",
        source_ref   => 'bowtie1_rRNAS_pm',
        seqcount_ref => [ "identical", ".dupcount\$" ],
        cqs_tools    => $def->{cqstools},
        option       => '',
        class        => 'CQS::BAMSequenceCountTable',
        prefix       => 'rRNAS_pm_'
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
        target_dir    => $def->{target_dir} . "/bowtie1_bacteria_group1_pm",
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
        target_dir   => $def->{target_dir} . "/bowtie1_bacteria_group1_pm_count",
        option       => $def->{smallrnacount_option},
        source_ref   => 'bowtie1_bacteria_group1_pm',
        cqs_tools    => $def->{cqstools},
        seqcount_ref => [ "identical", ".dupcount\$" ],
        'class'      => 'CQS::CQSChromosomeCount'
      },

      bowtie1_bacteria_group1_pm_table => {
        pbs => {
          'email'    => $def->{email},
          'walltime' => '72',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=1'
        },
        cluster    => $cluster,
        sh_direct  => 1,
        perform    => 1,
        target_dir => $def->{target_dir} . "/bowtie1_bacteria_group1_pm_table",
        source_ref => [ 'bowtie1_bacteria_group1_pm_count', '.xml' ],
        cqs_tools  => $def->{cqstools},
        option     => '',
        class      => 'CQS::CQSChromosomeTable',
        prefix     => 'bacteria_group1_pm_'
      },
      bowtie1_bacteria_group1_pm_table_vis => {
        class                => "CQS::UniqueR",
        perform              => 1,
        target_dir           => $def->{target_dir} . "/bowtie1_bacteria_group1_pm_table",
        rtemplate            => "bacteriaGroupMappingVis.R",
        output_file          => ".group1Mapping.Result",
        output_file_ext      => ".toSpecies.csv",
        parameterSampleFile1 => $groups,
        parameterFile1_ref   => [ "bowtie1_bacteria_group1_pm_table", ".count\$" ],
        parameterFile2       => $def->{bacteria_group1_log},
        sh_direct            => 1,
        pbs                  => {
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
        target_dir    => $def->{target_dir} . "/bowtie1_bacteria_group2_pm",
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
        target_dir   => $def->{target_dir} . "/bowtie1_bacteria_group2_pm_count",
        option       => $def->{smallrnacount_option},
        source_ref   => 'bowtie1_bacteria_group2_pm',
        cqs_tools    => $def->{cqstools},
        seqcount_ref => [ "identical", ".dupcount\$" ],
        'class'      => 'CQS::CQSChromosomeCount'
      },

      bowtie1_bacteria_group2_pm_table => {
        pbs => {
          'email'    => $def->{email},
          'walltime' => '72',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=1'
        },
        cluster    => $cluster,
        sh_direct  => 1,
        perform    => 1,
        target_dir => $def->{target_dir} . "/bowtie1_bacteria_group2_pm_table",
        source_ref => [ 'bowtie1_bacteria_group2_pm_count', '.xml' ],
        cqs_tools  => $def->{cqstools},
        option     => '',
        class      => 'CQS::CQSChromosomeTable',
        prefix     => 'bacteria_group2_pm_'
      },

      bowtie1_bacteria_group2_pm_table_vis => {
        class                => "CQS::UniqueR",
        perform              => 1,
        target_dir           => $def->{target_dir} . "/bowtie1_bacteria_group2_pm_table",
        rtemplate            => "bacteriaGroupMappingVis.R",
        output_file          => ".group2Mapping.Result",
        output_file_ext      => ".toSpecies.csv",
        parameterSampleFile1 => $groups,
        parameterFile1_ref   => [ "bowtie1_bacteria_group2_pm_table", ".count\$" ],
        parameterFile2       => $def->{bacteria_group2_log},
        sh_direct            => 1,
        pbs                  => {
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
        target_dir    => $def->{target_dir} . "/bowtie1_fungus_group4_pm",
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
        target_dir   => $def->{target_dir} . "/bowtie1_fungus_group4_pm_count",
        option       => $def->{smallrnacount_option},
        source_ref   => 'bowtie1_fungus_group4_pm',
        cqs_tools    => $def->{cqstools},
        seqcount_ref => [ "identical", ".dupcount\$" ],
        'class'      => 'CQS::CQSChromosomeCount'
      },

      bowtie1_fungus_group4_pm_table => {
        pbs => {
          'email'    => $def->{email},
          'walltime' => '72',
          'mem'      => '40gb',
          'nodes'    => '1:ppn=1'
        },
        cluster    => $cluster,
        sh_direct  => 1,
        perform    => 1,
        target_dir => $def->{target_dir} . "/bowtie1_fungus_group4_pm_table",
        source_ref => [ 'bowtie1_fungus_group4_pm_count', '.xml' ],
        cqs_tools  => $def->{cqstools},
        option     => '',
        class      => 'CQS::CQSChromosomeTable',
        prefix     => 'fungus_group4_pm_'
      },
      bowtie1_fungus_group4_pm_table_vis => {
        class                => "CQS::UniqueR",
        perform              => 1,
        target_dir           => $def->{target_dir} . "/bowtie1_fungus_group4_pm_table",
        rtemplate            => "bacteriaGroupMappingVis.R",
        output_file          => ".group4Mapping.Result",
        output_file_ext      => ".toSpecies.csv",
        parameterSampleFile1 => $groups,
        parameterFile1_ref   => [ "bowtie1_fungus_group4_pm_table", ".count\$" ],
        parameterFile2       => $def->{fungus_group4_log},
        sh_direct            => 1,
        pbs                  => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },
    };

    $config = merge( $config, $unmappedreads );

    push @individual,
      (
      "bowtie1_tRNA_pm",            "bowtie1_tRNA_pm_count",            "bowtie1_rRNAL_pm",           "bowtie1_rRNAL_pm_names",
      "bowtie1_rRNAS_pm",           "bowtie1_rRNAS_pm_names",           "bowtie1_bacteria_group1_pm", "bowtie1_bacteria_group1_pm_count",
      "bowtie1_bacteria_group2_pm", "bowtie1_bacteria_group2_pm_count", "bowtie1_fungus_group4_pm",   "bowtie1_fungus_group4_pm_count"
      );
    push @summary,
      (
      "bowtie1_tRNA_pm_table",            "bowtie1_tRNA_pm_table_vis",            "bowtie1_rRNAL_pm_table",           "bowtie1_rRNAS_pm_table",
      "bowtie1_bacteria_group1_pm_table", "bowtie1_bacteria_group1_pm_table_vis", "bowtie1_bacteria_group2_pm_table", "bowtie1_bacteria_group2_pm_table_vis",
      "bowtie1_fungus_group4_pm_table",   "bowtie1_fungus_group4_pm_table_vis",
      );

    push @mapped, ( "bowtie1_tRNA_pm_count", ".xml", "bowtie1_bacteria_group1_pm_count", ".xml", "bowtie1_bacteria_group2_pm_count", ".xml", "bowtie1_fungus_group4_pm_count", ".xml" );

    push @pmnames, ( "bowtie1_rRNAL_pm_names", "bowtie1_rRNAS_pm_names", );

    #do unmapped reads DESeq2
    if ($do_comparison) {
      my $unmapped_comparison = {

        #DESeq2
        group1_deseq2 => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $def->{target_dir} . "/bacteria_group1_deseq2",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_bacteria_group1_pm_table_vis", ".toSpecies.csv\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => 1,
          pvalue               => 0.05,
          fold_change          => 1.5,
          min_median_read      => 5,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
        group2_deseq2 => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $def->{target_dir} . "/bacteria_group2_deseq2",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_bacteria_group2_pm_table_vis", ".toSpecies.csv\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => 1,
          pvalue               => 0.05,
          fold_change          => 1.5,
          min_median_read      => 5,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
       group4_deseq2 => {
          class                => "Comparison::DESeq2",
          perform              => 1,
          target_dir           => $def->{target_dir} . "/fungus_group4_deseq2",
          option               => "",
          source_ref           => "pairs",
          groups_ref           => "groups",
          countfile_ref        => [ "bowtie1_fungus_group4_pm_table_vis", ".toSpecies.csv\$" ],
          sh_direct            => 1,
          show_DE_gene_cluster => 1,
          pvalue               => 0.05,
          fold_change          => 1.5,
          min_median_read      => 5,
          pbs                  => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "10",
            "mem"      => "10gb"
          },
        },
      };

      $config = merge( $config, $unmapped_comparison );
      push @summary, ( "group1_deseq2", "group2_deseq2" );
    }

  }

  if ($blast_unmapped_reads) {
    my $blast = {

      bowtie1_unmapped_reads => {
        class       => "CQS::Perl",
        perform     => 1,
        target_dir  => $def->{target_dir} . "/bowtie1_unmapped_reads",
        perlFile    => "unmappedReadsToFastq.pl",
        source_ref  => $identical_ref,
        source2_ref => \@mapped,
        source3_ref => \@pmnames,
        output_ext  => "_clipped_identical.unmapped.fastq.gz",
        output_other_ext  => "_clipped_identical.unmapped.fastq.dupcount",
        sh_direct   => 1,
        pbs         => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },

      bowtie1_unmapped_sequence_count_table => {
        class           => "CQS::SmallRNASequenceCountTable",
        perform         => 1,
        target_dir      => $def->{target_dir} . "/bowtie1_unmapped_sequence_count_table",
        option          => "",
        source_ref      => [ "identical", ".dupcount\$" ],
        fastq_files_ref => [ "bowtie1_unmapped_reads", ".fastq.gz" ],
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
        target_dir => $def->{target_dir} . "/bowtie1_unmapped_sequence_blast",
        option     => "",
        source_ref => [ "bowtie1_unmapped_sequence_count_table", ".fasta\$" ],
        sh_direct  => 0,
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
    push @individual, ("bowtie1_unmapped_reads");
    push @summary, ( "bowtie1_unmapped_sequence_count_table", "bowtie1_unmapped_sequence_blast" );
  }
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
      "walltime" => "72",
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

  my $config = getParclipSmallRNAConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
