#!/usr/bin/perl
package Pipeline::SmallRNA;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
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
  my $task_name = $def->{task_name};

  my $search_not_identical   = getValue( $def, "search_not_identical" );
  my $search_host_genome     = getValue( $def, "search_host_genome" );
  my $search_nonhost_genome  = getValue( $def, "search_nonhost_genome" );
  my $search_nonhost_library = getValue( $def, "search_nonhost_library" );
  my $search_nonhost_database = $search_nonhost_genome || $search_nonhost_library;

  my $blast_top_reads      = getValue( $def, "blast_top_reads" );
  my $blast_unmapped_reads = getValue( $def, "blast_unmapped_reads" );

  my $top_read_number = getValue( $def, "top_read_number" );
  my $host_genome_dir;
  if ($search_host_genome) {
    $host_genome_dir = create_directory_or_die( $def->{target_dir} . "/host_genome" );
  }

  my $nonhost_library_dir;
  if ($search_nonhost_library) {
    $nonhost_library_dir = create_directory_or_die( $def->{target_dir} . "/nonhost_library" );
  }

  my $nonhost_genome_dir;
  if ($search_nonhost_genome) {
    $nonhost_genome_dir = create_directory_or_die( $def->{target_dir} . "/nonhost_genome" );
  }

  my $nonhost_blast_dir;
  if ( $blast_unmapped_reads || $search_nonhost_database ) {
    $nonhost_blast_dir = create_directory_or_die( $def->{target_dir} . "/final_unmapped" );
  }

  my $data_visualization_dir = create_directory_or_die( $def->{target_dir} . "/data_visualization" );

  my @table_for_correlation = ( "identical_sequence_count_table", "^(?!.*?read).*\.count\$" );
  my @table_for_countSum    = ();
  my @table_for_pieSummary  = ();
  my @name_for_pieSummary   = ();
  my @table_for_readSummary = ();
  my @name_for_readSummary  = ();

  #print Dumper($config);

  my $do_comparison     = defined $def->{pairs};
  my $groups            = $def->{groups};
  my $groups_vis_layout = $def->{groups_vis_layout};
  if ( !defined $def->{groups_smallRNA_vis_layout} ) {
    $def->{groups_smallRNA_vis_layout} = $def->{groups_vis_layout};
  }

  my $DE_min_median_read_top      = getValue( $def, "DE_min_median_read_top" );
  my $DE_min_median_read_smallRNA = getValue( $def, "DE_min_median_read_smallRNA" );

  my $max_sequence_extension_base = getValue( $def, "max_sequence_extension_base" );
  $def->{nonhost_table_option} = "--maxExtensionBase " . $def->{max_sequence_extension_base} . " " . $def->{nonhost_table_option};
  my $perform_contig_analysis = $def->{perform_contig_analysis};
  if ($perform_contig_analysis) {
    $def->{nonhost_table_option} = $def->{nonhost_table_option} . " --outputReadContigTable";
  }

  my $deseq2Task;
  my $bowtie1Task;
  my $bowtie1CountTask;
  my $bowtie1TableTask;

  if ($do_comparison) {
    my $taskKey = "top${top_read_number}";

    $deseq2Task = addDEseq2( $config, $def, $summary_ref, "${taskKey}_reads", [ "identical_sequence_count_table", ".read.count\$" ], $class_independent_dir, $DE_min_median_read_top );
    addDeseq2Visualization( $config, $def, $summary_ref, "${taskKey}_reads", [ $deseq2Task, "_DESeq2.csv\$" ], $data_visualization_dir, "pairs_top_deseq2_vis_layout" );

    $deseq2Task = addDEseq2( $config, $def, $summary_ref, "${taskKey}_contigs", [ "identical_sequence_count_table", ".count\$" ], $class_independent_dir, $DE_min_median_read_top );
    addDeseq2Visualization( $config, $def, $summary_ref, "${taskKey}_contigs", [ $deseq2Task, "_DESeq2.csv\$" ], $data_visualization_dir, "pairs_top_deseq2_vis_layout" );

    $deseq2Task = addDEseq2( $config, $def, $summary_ref, "${taskKey}_minicontigs", [ "identical_sequence_count_table", ".minicontig.count\$" ], $class_independent_dir, $DE_min_median_read_top );
    addDeseq2Visualization( $config, $def, $summary_ref, "${taskKey}_minicontigs", [ $deseq2Task, "_DESeq2.csv\$" ], $data_visualization_dir, "pairs_top_deseq2_vis_layout" );
  }

  my $identical_ref = [ "identical", ".fastq.gz\$" ];

  if ($search_host_genome) {
    getValue( $def, "coordinate" );

    #1 mismatch search, NTA
    addBowtie( $config, $def, $individual_ref, "bowtie1_genome_1mm_NTA", $host_genome_dir, $def->{bowtie1_index}, [ "identical_NTA", ".fastq.gz\$" ], $def->{bowtie1_option_1mm} );

    my $host_genome = {

      bowtie1_genome_1mm_NTA_smallRNA_count => {
        class           => "CQS::SmallRNACount",
        perform         => 1,
        target_dir      => $host_genome_dir . "/bowtie1_genome_1mm_NTA_smallRNA_count",
        option          => $def->{host_smallrnacount_option},
        source_ref      => "bowtie1_genome_1mm_NTA",
        fastq_files_ref => "identical_NTA",
        seqcount_ref    => [ "identical", ".dupcount\$" ],
        cqs_tools       => $def->{cqstools},
        coordinate_file => $def->{coordinate},
        fasta_file      => $def->{coordinate_fasta},
        sh_direct       => 1,
        cluster         => $cluster,
        pbs             => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "72",
          "mem"       => "40gb"
        },
      },
      bowtie1_genome_1mm_NTA_smallRNA_table => {
        class      => "CQS::SmallRNATable",
        perform    => 1,
        target_dir => $host_genome_dir . "/bowtie1_genome_1mm_NTA_smallRNA_table",
        option     => $def->{host_smallrnacounttable_option},
        source_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".mapped.xml" ],
        cqs_tools  => $def->{cqstools},
        prefix     => "smallRNA_1mm_",
        hasYRNA    => $def->{hasYRNA},
        sh_direct  => 1,
        cluster    => $cluster,
        pbs        => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "10",
          "mem"       => "10gb"
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
        parameterSampleFile3      => $def->{groups_smallRNA_vis_layout},
        rCode                     => 'textSize=9;groupTextSize=' . $def->{table_vis_group_text_size} . ';',
        sh_direct                 => 1,
        pbs                       => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
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
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      },
    };

    if ( $def->{has_NTA} && $def->{consider_tRNA_NTA} ) {
      $host_genome->{"bowtie1_genome_1mm_NTA_smallRNA_count"}{"cca_file_ref"} = "identical_check_cca";
    }

    push @table_for_pieSummary, ( "bowtie1_genome_1mm_NTA_smallRNA_count", ".count\$" );
    push @name_for_pieSummary, "Host Small RNA";
    push @table_for_correlation, ( "bowtie1_genome_1mm_NTA_smallRNA_table", "^(?!.*?read).*\.count\$" );
    push @table_for_readSummary,
      ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.read.count\$", "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.read.count\$", "bowtie1_genome_1mm_NTA_smallRNA_table", ".other.read.count\$" );
    push @name_for_readSummary, ( "Host miRNA", "Host tRNA", "Host other small RNA" );
    push @table_for_countSum,
      ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.count\$", "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.count\$", "bowtie1_genome_1mm_NTA_smallRNA_table", ".other.count\$" );
    push @$individual_ref, ( "bowtie1_genome_1mm_NTA", "bowtie1_genome_1mm_NTA_smallRNA_count" );
    push @$summary_ref, ( "bowtie1_genome_1mm_NTA_smallRNA_table", "bowtie1_genome_1mm_NTA_smallRNA_category", "host_genome_tRNA_category" );

    $config = merge( $config, $host_genome );
    if ($do_comparison) {
      my @visual_source = ();

      #miRNA
      $deseq2Task = addDEseq2( $config, $def, \@$summary_ref, "miRNA", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.count\$" ], $host_genome_dir, $DE_min_median_read_smallRNA );
      push( @visual_source, ( $deseq2Task, "_DESeq2.csv\$" ) );
      addDEseq2( $config, $def, $summary_ref, "miRNA_NTA",        [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.NTA.count\$" ],        $host_genome_dir, $DE_min_median_read_smallRNA );
      addDEseq2( $config, $def, $summary_ref, "miRNA_NTA_base",   [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.NTA.base.count\$" ],   $host_genome_dir, $DE_min_median_read_smallRNA );
      addDEseq2( $config, $def, $summary_ref, "miRNA_isomiR",     [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.isomiR.count\$" ],     $host_genome_dir, $DE_min_median_read_smallRNA );
      addDEseq2( $config, $def, $summary_ref, "miRNA_isomiR_NTA", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.isomiR_NTA.count\$" ], $host_genome_dir, $DE_min_median_read_smallRNA );
      addDEseq2( $config, $def, $summary_ref, "miRNA_reads",      [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.read.count\$" ],       $host_genome_dir, $DE_min_median_read_smallRNA );
      addDeseq2Visualization( $config, $def, $summary_ref, "host_genome_miRNA",
        [ "deseq2_miRNA_isomiR", "_DESeq2.csv\$", "deseq2_miRNA_NTA", "_DESeq2.csv\$", "deseq2_miRNA_isomiR_NTA", "_DESeq2.csv\$" ],
        $data_visualization_dir, "pairs_host_miRNA_deseq2_vis_layout" );

      #tRNA
      $deseq2Task = addDEseq2( $config, $def, $summary_ref, "tRNA", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.count\$" ], $host_genome_dir, $DE_min_median_read_smallRNA );
      push( @visual_source, ( $deseq2Task, "_DESeq2.csv\$" ) );
      addDEseq2( $config, $def, $summary_ref, "tRNA_reads",     [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.read.count\$" ],      $host_genome_dir, $DE_min_median_read_smallRNA );
      addDEseq2( $config, $def, $summary_ref, "tRNA_aminoacid", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.aminoacid.count\$" ], $host_genome_dir, $DE_min_median_read_smallRNA );

      if ( $def->{host_smallrnacounttable_option} =~ /yRNAsnRNAsnoRNA/ ) {
        if ( $def->{hasYRNA} ) {

          #yRNA
          $deseq2Task = addDEseq2( $config, $def, $summary_ref, "yRNA", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".yRNA.count\$" ], $host_genome_dir, $DE_min_median_read_smallRNA );
          push( @visual_source, ( $deseq2Task, "_DESeq2.csv\$" ) );
          addDEseq2( $config, $def, $summary_ref, "yRNA_reads", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".yRNA.read.count\$" ], $host_genome_dir, $DE_min_median_read_smallRNA );
        }

        #snRNA
        $deseq2Task = addDEseq2( $config, $def, $summary_ref, "snRNA", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".snRNA.count\$" ], $host_genome_dir, $DE_min_median_read_smallRNA );
        push( @visual_source, ( $deseq2Task, "_DESeq2.csv\$" ) );
        addDEseq2( $config, $def, $summary_ref, "snRNA_reads", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".snRNA.read.count\$" ], $host_genome_dir, $DE_min_median_read_smallRNA );

        #snoRNA
        $deseq2Task = addDEseq2( $config, $def, $summary_ref, "snoRNA", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".snoRNA.count\$" ], $host_genome_dir, $DE_min_median_read_smallRNA );
        push( @visual_source, ( $deseq2Task, "_DESeq2.csv\$" ) );
        addDEseq2( $config, $def, $summary_ref, "snoRNA_reads", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".snoRNA.read.count\$" ], $host_genome_dir, $DE_min_median_read_smallRNA );
      }

      #otherSmallRNA
      $deseq2Task = addDEseq2( $config, $def, $summary_ref, "otherSmallRNA", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".other.count\$" ], $host_genome_dir, $DE_min_median_read_smallRNA );
      push( @visual_source, ( $deseq2Task, "_DESeq2.csv\$" ) );
      addDEseq2( $config, $def, $summary_ref, "otherSmallRNA_reads", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".other.read.count\$" ], $host_genome_dir, $DE_min_median_read_smallRNA );

      #host genome smallRNA visualization
      addDeseq2Visualization( $config, $def, $summary_ref, "host_genome", \@visual_source, $data_visualization_dir, "pairs_host_deseq2_vis_layout" );
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
        rtemplate                => "countTableVisFunctions.R,tRNAPositionVis.R",
        output_file              => ".tRNAPositionVis",
        output_file_ext          => ".alltRNAPosition.png",
        parameterSampleFile1_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".tRNA.position\$" ],
        parameterSampleFile2     => $trna_vis_groups,
        parameterSampleFile3_ref => $trna_sig_result,
        sh_direct                => 1,
        pbs                      => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push @$summary_ref, ("host_genome_tRNA_PositionVis");
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
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "10",
          "mem"       => "10gb"
        },
      },

      bowtie1_genome_unmapped_reads => {
        class       => "CQS::Perl",
        perform     => 1,
        target_dir  => $host_genome_dir . "/bowtie1_genome_unmapped_reads",
        perlFile    => "unmappedReadsToFastq.pl",
        source_ref  => [ "identical", ".fastq.gz\$" ],
        source2_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".mapped.xml" ],
        source3_ref => ["bowtie1_genome_1mm_NTA_pmnames"],
        output_ext  => "_clipped_identical.unmapped.fastq.gz",
        output_other_ext =>
"_clipped_identical.unmapped.fastq.dupcount,_clipped_identical.mappedToHostGenome.dupcount,_clipped_identical.mappedToHostGenome.fastq.gz,_clipped_identical.short.fastq.gz,_clipped_identical.short.dupcount",
        sh_direct => 1,
        pbs       => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      },
      bowtie1_genome_host_reads_table => {
        class      => "CQS::CQSDatatable",
        perform    => 1,
        target_dir => $host_genome_dir . "/bowtie1_genome_host_reads_table",
        source_ref => [ "bowtie1_genome_unmapped_reads", ".mappedToHostGenome.dupcount\$" ],
        option     => "-k 2 -v 1 --fillMissingWithZero",
        cqstools   => $def->{cqstools},
        sh_direct  => 1,
        pbs        => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      }
    };
    $config = merge( $config, $unmapped_reads );
    push @$individual_ref, ( "bowtie1_genome_1mm_NTA_pmnames", "bowtie1_genome_unmapped_reads" );
    push @$summary_ref, ("bowtie1_genome_host_reads_table");
    push @table_for_pieSummary,
      ( "bowtie1_genome_unmapped_reads", ".mappedToHostGenome.dupcount", "bowtie1_genome_unmapped_reads", ".short.dupcount", "bowtie1_genome_unmapped_reads", ".unmapped.fastq.dupcount" );
    push @name_for_pieSummary, ( "Mapped to Host Genome", "Too Short for Mapping", "Unmapped In Host" );
    push @table_for_readSummary, ( "bowtie1_genome_host_reads_table", ".count\$" );
    push @name_for_readSummary, ("Host Genome");
    $identical_ref = [ "bowtie1_genome_unmapped_reads", ".unmapped.fastq.gz\$" ];
  }

  my @mapped  = ();
  my @pmnames = ();
  my @overlap = ();

  if ($search_nonhost_database) {

    #Mapping host genome reads to non-host databases
    addNonhostDatabase(
      $config, $def, $individual_ref, $summary_ref, "HostGenomeReads_NonHost_pm", $nonhost_library_dir,    #general option
      $def->{bowtie1_all_nonHost_index}, [ "bowtie1_genome_unmapped_reads", ".mappedToHostGenome.fastq.gz" ],    #bowtie option
      $def->{smallrnacount_option} . ' --keepChrInName --categoryMapFile ' . $def->{all_nonHost_map},            #count option
      $def->{nonhost_table_option}                                                                               #table option
    );
    addNonhostVis(
      $config, $def,
      $summary_ref,
      "HostGenomeReads_NonHost_vis",
      $data_visualization_dir,
      {
        rtemplate          => "countTableVisFunctions.R,countTableVis.R",
        output_file        => ".NonHostAll.Result",
        output_file_ext    => ".Barplot.png",
        parameterFile1_ref => [ "bowtie1_HostGenomeReads_NonHost_pm_table", ".count\$" ],
      }
    );
  }

  #Mapping unmapped reads to nonhost genome
  if ($search_nonhost_genome) {
    addNonhostDatabase(
      $config, $def, $individual_ref, $summary_ref, "bacteria_group1_pm", $nonhost_genome_dir,    #general option
      $def->{bowtie1_bacteria_group1_index}, $identical_ref,                                      #bowtie option
      $def->{smallrnacount_option} . ' --keepChrInName --keepSequence',                           #count option
      $def->{nonhost_table_option} . ' --categoryMapFile ' . $def->{bacteria_group1_species_map}  #table option
    );
    addNonhostVis(
      $config, $def,
      $summary_ref,
      "nonhost_genome_bacteria_group1_vis",
      $data_visualization_dir,
      {
        rtemplate          => "countTableVisFunctions.R,countTableVis.R",
        output_file        => ".group1Mapping.Result",
        output_file_ext    => ".Piechart.png",
        parameterFile1_ref => [ "bowtie1_bacteria_group1_pm_table", ".category.count\$" ],
      }
    );

    #Mapping unmapped reads to group2 database
    addNonhostDatabase(
      $config, $def, $individual_ref, $summary_ref, "bacteria_group2_pm", $nonhost_genome_dir,    #general option
      $def->{bowtie1_bacteria_group2_index}, $identical_ref,                                      #bowtie option
      $def->{smallrnacount_option} . ' --keepChrInName --keepSequence',                           #count option
      $def->{nonhost_table_option} . ' --categoryMapFile ' . $def->{bacteria_group2_species_map}  #table option
    );
    addNonhostVis(
      $config, $def,
      $summary_ref,
      "nonhost_genome_bacteria_group2_vis",
      $data_visualization_dir,
      {
        rtemplate          => "countTableVisFunctions.R,countTableVis.R",
        output_file        => ".group2Mapping.Result",
        output_file_ext    => ".Piechart.png",
        parameterFile1_ref => [ "bowtie1_bacteria_group2_pm_table", ".category.count\$" ],
      }
    );

    #Mapping unmapped reads to group4 database
    addNonhostDatabase(
      $config, $def, $individual_ref, $summary_ref, "fungus_group4_pm", $nonhost_genome_dir,    #general option
      $def->{bowtie1_fungus_group4_index}, $identical_ref,                                      #bowtie option
      $def->{smallrnacount_option} . ' --keepChrInName --keepSequence',                         #count option
      $def->{nonhost_table_option} . ' --categoryMapFile ' . $def->{fungus_group4_species_map}  #table option
    );
    addNonhostVis(
      $config, $def,
      $summary_ref,
      "nonhost_genome_fungus_group4_vis",
      $data_visualization_dir,
      {
        rtemplate          => "countTableVisFunctions.R,countTableVis.R",
        output_file        => ".group4Mapping.Result",
        output_file_ext    => ".Piechart.png",
        parameterFile1_ref => [ "bowtie1_fungus_group4_pm_table", ".category.count\$" ],
      }
    );

    push @table_for_correlation,
      ( "bowtie1_bacteria_group1_pm_table", ".category.count\$", "bowtie1_bacteria_group2_pm_table", ".category.count\$", "bowtie1_fungus_group4_pm_table", ".category.count\$" );
    push @table_for_countSum,
      ( "bowtie1_bacteria_group1_pm_table", ".category.count\$", "bowtie1_bacteria_group2_pm_table", ".category.count\$", "bowtie1_fungus_group4_pm_table", ".category.count\$" );
    push @table_for_readSummary, ( "bowtie1_bacteria_group1_pm_table", ".read.count\$", "bowtie1_bacteria_group2_pm_table", ".read.count\$", "bowtie1_fungus_group4_pm_table", ".read.count\$" );
    push @name_for_readSummary, ( "Human Microbiome Bacteria", "Environment Bacteria", "Fungus" );

    push @mapped,  ( "bowtie1_bacteria_group1_pm_count", ".xml",          "bowtie1_bacteria_group2_pm_count", ".xml",          "bowtie1_fungus_group4_pm_count", ".xml" );
    push @overlap, ( "bowtie1_bacteria_group1_pm_table", ".read.count\$", "bowtie1_bacteria_group2_pm_table", ".read.count\$", "bowtie1_fungus_group4_pm_table", ".read.count\$" );

    if ($do_comparison) {

      addDEseq2( $config, $def, $summary_ref, "bacteria_group1",       [ "bowtie1_bacteria_group1_pm_table", ".category.count\$" ], $nonhost_genome_dir, $DE_min_median_read_smallRNA );
      addDEseq2( $config, $def, $summary_ref, "bacteria_group1_reads", [ "bowtie1_bacteria_group1_pm_table", ".read.count\$" ],     $nonhost_genome_dir, $DE_min_median_read_smallRNA );
      addDEseq2( $config, $def, $summary_ref, "bacteria_group2",       [ "bowtie1_bacteria_group2_pm_table", ".category.count\$" ], $nonhost_genome_dir, $DE_min_median_read_smallRNA );
      addDEseq2( $config, $def, $summary_ref, "bacteria_group2_reads", [ "bowtie1_bacteria_group2_pm_table", ".read.count\$" ],     $nonhost_genome_dir, $DE_min_median_read_smallRNA );
      addDEseq2( $config, $def, $summary_ref, "fungus_group4",         [ "bowtie1_fungus_group4_pm_table",   ".category.count\$" ], $nonhost_genome_dir, $DE_min_median_read_smallRNA );
      addDEseq2( $config, $def, $summary_ref, "fungus_group4_reads",   [ "bowtie1_fungus_group4_pm_table",   ".read.count\$" ],     $nonhost_genome_dir, $DE_min_median_read_smallRNA );

      addDeseq2Visualization( $config, $def, $summary_ref, "nonhost_genome_deseq2",
        [ "deseq2_bacteria_group1", "_DESeq2.csv\$", "deseq2_bacteria_group2", "_DESeq2.csv\$", "deseq2_fungus_group4", "_DESeq2.csv\$" ],
        $data_visualization_dir, "pairs_nonHostGroups_deseq2_vis_layout" );

      addDeseq2Visualization( $config, $def, $summary_ref, "nonhost_genome_deseq2_reads",
        [ "deseq2_bacteria_group1_reads", "_DESeq2.csv\$", "deseq2_bacteria_group2_reads", "_DESeq2.csv\$", "deseq2_fungus_group4_reads", "_DESeq2.csv\$" ],
        $data_visualization_dir, "pairs_nonHostGroups_deseq2_vis_layout" );
    }
  }

  #Mapping unmapped reads to nonhost library
  if ($search_nonhost_library) {

    #Mapping unmapped reads to miRBase library
    addNonhostDatabase(
      $config, $def, $individual_ref, $summary_ref, "miRBase_pm", $nonhost_library_dir,    #general option
      $def->{bowtie1_miRBase_index}, $identical_ref,                                       #bowtie option
      $def->{mirbase_count_option} . " -m --keepChrInName --keepSequence",                 #count option
      $def->{nonhost_table_option}                                                         #table option
    );

    push @table_for_countSum, ( "bowtie1_miRBase_pm_table", "^(?!.*?read).*\.count\$" );
    push @mapped,             ( "bowtie1_miRBase_pm_count", ".xml" );

    #Mapping unmapped reads to tRNA library
    addNonhostDatabase(
      $config, $def, $individual_ref, $summary_ref, "tRNA_pm", $nonhost_library_dir,       #general option
      $def->{bowtie1_tRNA_index}, $identical_ref,                                          #bowtie option
      $def->{smallrnacount_option} . " --keepChrInName --keepSequence",                    #count option
      $def->{nonhost_table_option} . ' --categoryMapFile ' . $def->{trna_category_map}     #table option
    );
    addNonhostVis(
      $config, $def,
      $summary_ref,
      "nonhost_library_tRNA_vis",
      $data_visualization_dir,
      {
        rtemplate          => "countTableVisFunctions.R,bacteriaTrnaMappingVis.R",
        output_file        => ".tRNAMapping.Result",
        output_file_ext    => ".Species12.csv;.tRNAType1.csv;.tRNAType2.csv",
        parameterFile1_ref => [ "bowtie1_tRNA_pm_table", ".count\$" ],
      }
    );

    #Mapping unmapped reads to rRNA library
    addNonhostDatabase(
      $config, $def, $individual_ref, $summary_ref, "rRNA_pm", $nonhost_library_dir,    #general option
      $def->{bowtie1_rRNA_index}, $identical_ref,                                       #bowtie option
      $def->{smallrnacount_option} . ' --keepChrInName --keepSequence --categoryMapFile ' . $def->{rrna_category_map},    #count option                                          #count option
      $def->{nonhost_table_option}                                                                                        #table option
    );
    addNonhostVis(
      $config, $def,
      $summary_ref,
      "nonhost_library_rRNA_vis",
      $data_visualization_dir,
      {
        rtemplate          => "countTableVisFunctions.R,countTableVis.R",
        output_file        => ".rRNAMapping.Result",
        output_file_ext    => ".Barplot.png",
        parameterFile1_ref => [ "bowtie1_rRNA_pm_table", ".count\$" ],
      }
    );

    push @table_for_correlation, ( "bowtie1_tRNA_pm_table", "^(?!.*?read).*\.count\$", "bowtie1_rRNA_pm_table", "^(?!.*?read).*\.count\$" );
    push @table_for_countSum,    ( "bowtie1_tRNA_pm_table", ".category.count\$",       "bowtie1_rRNA_pm_table", "$task_name\.count\$" );
    push @table_for_readSummary, ( "bowtie1_tRNA_pm_table", ".read.count\$",           "bowtie1_rRNA_pm_table", ".read.count\$" );
    push @name_for_readSummary,  ( "Non host tRNA",         "Non host rRNA" );
    push @mapped,                ( "bowtie1_tRNA_pm_count", ".xml",                    "bowtie1_rRNA_pm_count", ".xml" );
    push @overlap,               ( "bowtie1_tRNA_pm_table", ".read.count\$",           "bowtie1_rRNA_pm_table", ".read.count\$" );

    if ($do_comparison) {

      addDEseq2( $config, $def, $summary_ref, "nonhost_tRNA",           [ "bowtie1_tRNA_pm_table",    ".count\$" ],          $nonhost_library_dir, $DE_min_median_read_smallRNA );
      addDEseq2( $config, $def, $summary_ref, "nonhost_tRNA_reads",     [ "bowtie1_tRNA_pm_table",    ".read.count\$" ],     $nonhost_library_dir, $DE_min_median_read_smallRNA );
      addDEseq2( $config, $def, $summary_ref, "nonhost_tRNA_category",  [ "bowtie1_tRNA_pm_table",    ".category.count\$" ], $nonhost_library_dir, $DE_min_median_read_smallRNA );
      addDEseq2( $config, $def, $summary_ref, "nonhost_tRNA_species",   [ "nonhost_library_tRNA_vis", ".Species12.csv\$" ],  $nonhost_library_dir, $DE_min_median_read_smallRNA );
      addDEseq2( $config, $def, $summary_ref, "nonhost_tRNA_type",      [ "nonhost_library_tRNA_vis", ".tRNAType1.csv\$" ],  $nonhost_library_dir, $DE_min_median_read_smallRNA );
      addDEseq2( $config, $def, $summary_ref, "nonhost_tRNA_anticodon", [ "nonhost_library_tRNA_vis", ".tRNAType2.csv\$" ],  $nonhost_library_dir, $DE_min_median_read_smallRNA );

      addDeseq2Visualization(
        $config, $def,
        $summary_ref,
        "nonhost_library_deseq2",
        [
          "deseq2_nonhost_tRNA",      "_DESeq2.csv\$", "deseq2_nonhost_tRNA_category",  "_DESeq2.csv\$", "deseq2_nonhost_tRNA_species", "_DESeq2.csv\$",
          "deseq2_nonhost_tRNA_type", "_DESeq2.csv\$", "deseq2_nonhost_tRNA_anticodon", "_DESeq2.csv\$"
        ],
        $data_visualization_dir,
        "pairs_nonHostLibrary_deseq2_vis_layout"
      );

      addDEseq2( $config, $def, $summary_ref, "nonhost_rRNA", [ "bowtie1_rRNA_pm_table", ".count\$" ], $nonhost_library_dir, $DE_min_median_read_smallRNA );
    }
  }

  if ($search_nonhost_database) {

    $config->{nonhost_overlap_vis} = {
      class                     => "CQS::UniqueR",
      perform                   => 1,
      target_dir                => $data_visualization_dir . "/nonhost_overlap_vis",
      rtemplate                 => "countTableVisFunctions.R,NonHostOverlap.R",
      output_file               => ".NonHost.Reads",
      output_file_ext           => ".Overlap.csv",
      parameterSampleFile1_ref  => \@overlap,
      parameterSampleFile2Order => $def->{groups_order},
      parameterSampleFile2      => $groups,
      parameterSampleFile3      => $groups_vis_layout,
      parameterFile3_ref        => [ "fastqc_count_vis", ".Reads.csv\$" ],
      sh_direct                 => 1,
      rCode                     => 'maxCategory=8;textSize=9;groupTextSize=' . $def->{table_vis_group_text_size} . ';',
      pbs                       => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    };
    push @$summary_ref, "nonhost_overlap_vis";

    $config->{final_unmapped_reads} = {
      class            => "CQS::Perl",
      perform          => 1,
      target_dir       => $nonhost_blast_dir . "/final_unmapped_reads",
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
    };

    $identical_ref = [ "final_unmapped_reads", ".fastq.gz\$" ];
    push @$individual_ref,      ("final_unmapped_reads");
    push @table_for_pieSummary, ( "final_unmapped_reads", ".dupcount" );
    push @name_for_pieSummary,  ("UnMapped");
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
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
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
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "12",
      "mem"       => "10gb"
    },
    };

  my $name_for_pieSummary_r = "readFilesModule=c('" . join( "','", @name_for_pieSummary ) . "')";
  $config->{reads_in_tasks_pie} = {
    class                    => "CQS::UniqueR",
    suffix                   => "_pie",
    perform                  => 1,
    target_dir               => $data_visualization_dir . "/reads_in_tasks",
    rtemplate                => "countTableVisFunctions.R,ReadsInTasksPie.R",
    output_file_ext          => ".NonParallel.TaskReads.csv",
    parameterSampleFile1_ref => \@table_for_pieSummary,
    parameterSampleFile2     => $groups,
    parameterSampleFile3     => $groups_vis_layout,
    rCode                    => $name_for_pieSummary_r,

    #    parameterFile3_ref       => [ "fastqc_count_vis", ".Reads.csv\$" ],
    sh_direct => 1,
    pbs       => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "12",
      "mem"       => "10gb"
    },
  };
  my $name_for_readSummary_r = "readFilesModule=c('" . join( "','", @name_for_readSummary ) . "')";
  $config->{reads_mapping_summary} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $data_visualization_dir . "/reads_mapping_summary",
    rtemplate                => "countTableVisFunctions.R,ReadsMappingSummary.R",
    output_file_ext          => ".ReadsMapping.Summary.csv",
    parameterFile1_ref       => [ "identical_sequence_count_table", $task_name . "_sequence.read.count\$" ],
    parameterSampleFile1_ref => \@table_for_readSummary,
    parameterSampleFile2     => $groups,
    parameterSampleFile3     => $groups_vis_layout,
    rCode                    => $name_for_readSummary_r,
    sh_direct                => 1,
    pbs                      => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "12",
      "mem"       => "10gb"
    },
  };
  push @$summary_ref, ( "count_table_correlation", "reads_in_tasks", "reads_in_tasks_pie", "reads_mapping_summary" );

  #add time cost task in the end of pipeline
  #search not identical reads to genome, for IGV
  if ( $search_host_genome && $search_not_identical ) {
    addBowtie( $config, $def, $individual_ref, "bowtie1_genome_1mm_notidentical", $host_genome_dir, $def->{bowtie1_index}, $not_identical_ref, $def->{bowtie1_option_1mm} );
  }

  #blast top reads
  if ($blast_top_reads) {
    if ($do_comparison) {
      addDeseq2SignificantSequenceBlastn( $config, $def, $summary_ref, "deseq2_top${top_read_number}_minicontigs", $class_independent_dir );
    }
    else {
      #addBlastn( $config, $def, $summary_ref, "identical_sequence_top${top_read_number}_contig_blast",     "identical_sequence_count_table", "sequence.count.fasta\$",   $class_independent_dir );
      #addBlastn( $config, $def, $summary_ref, "identical_sequence_top${top_read_number}_read_blast",       "identical_sequence_count_table", "read.count.fasta\$",       $class_independent_dir );
      addBlastn( $config, $def, $summary_ref, "identical_sequence_top${top_read_number}_minicontig_blast", "identical_sequence_count_table", "minicontig.count.fasta\$", $class_independent_dir );
    }
  }

  #blast unmapped reads
  if ($blast_unmapped_reads) {
    $config->{"unmapped_sequence_count_table"} = {
      class      => "CQS::SmallRNASequenceCountTable",
      perform    => 1,
      target_dir => $nonhost_blast_dir . "/unmapped_sequence_count_table",
      option     => "--maxExtensionBase $max_sequence_extension_base -n $top_read_number --exportFastaNumber $top_read_number",

      source_ref      => [ "identical", ".dupcount\$" ],
      fastq_files_ref => $identical_ref,
      cqs_tools       => $def->{cqstools},
      suffix          => "_unmapped",
      sh_direct       => 1,
      cluster         => $cluster,
      pbs             => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "10gb"
      },
    };
    push @$summary_ref, "unmapped_sequence_count_table";

    if ($do_comparison) {
      $deseq2Task = addDEseq2( $config, $def, $summary_ref, "unmapped_sequence_minicontigs", [ "unmapped_sequence_count_table", ".minicontig.count\$" ], $nonhost_blast_dir, $DE_min_median_read_top );
      addDeseq2SignificantSequenceBlastn( $config, $def, $summary_ref, $deseq2Task, $nonhost_blast_dir );
    }
    else {
      addBlastn( $config, $def, $summary_ref, "unmapped_sequence_minicontig_blast", "unmapped_sequence_count_table", "minicontig.count.fasta\$", $nonhost_blast_dir );
    }
  }

  $config->{sequencetask} = {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => $def->{target_dir} . "/sequencetask",
    option     => "",
    source     => {
      step1 => $individual_ref,
      step2 => $summary_ref,
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
