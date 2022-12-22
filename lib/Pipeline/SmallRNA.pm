#!/usr/bin/perl
package Pipeline::SmallRNA;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use CQS::StringUtils;
use Pipeline::PipelineUtils;
use Pipeline::SmallRNAUtils;
use Data::Dumper;
use Hash::Merge qw( merge );
use Storable qw(dclone);

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performSmallRNA performSmallRNATask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.06';

sub getUniqueGroups {
  my ($def) = @_;
  my $result = "";
  if (defined $def->{unique_groups}){
    my $uni_groups = $def->{unique_groups};
    my @uniqueGroups;
    if (is_string($uni_groups)){
      @uniqueGroups = split(/[,;]/, $uni_groups);
    }else{
      @uniqueGroups = @$uni_groups;
    }
    $result = "uniqueGroupNames=c('" . join("','", @uniqueGroups) . "');";
  }
  return($result);
}

sub getSmallRNAConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  initializeSmallRNADefaultOptions($def);

  my ( $config, $individual_ref, $summary_ref, $cluster, $not_identical_ref, $preprocessing_dir, $class_independent_dir, $identical_ref, $host_identical_ref ) = getPrepareConfig( $def, 1 );

  #merge summary and individual 
  push @$individual_ref, @$summary_ref;
  $summary_ref = $individual_ref;

  my $task_name = $def->{task_name};

  my $hasMicroRNAOnly = getValue( $def, "hasMicroRNAOnly", 0 );
  my $notMicroRNAOnly = !$hasMicroRNAOnly;

  my $search_not_identical   = getValue( $def, "search_not_identical" );
  my $search_host_genome     = getValue( $def, "search_host_genome" );
  my $search_nonhost_genome  = getValue( $def, "search_nonhost_genome" ) && $notMicroRNAOnly;
  my $search_nonhost_library = getValue( $def, "search_nonhost_library" ) && $notMicroRNAOnly;

  my $search_refseq_genome     = getValue( $def, "search_refseq_genome", 0 ) && $notMicroRNAOnly;

  my $search_nonhost_database = $search_nonhost_genome || $search_nonhost_library || $search_refseq_genome;

  my $perform_annotate_unmapped_reads    = getValue( $def, "perform_annotate_unmapped_reads" );
  my $perform_class_independent_analysis = getValue( $def, "perform_class_independent_analysis", 1 );
  my $perform_short_reads_source         = $search_host_genome && getValue( $def, "perform_short_reads_source" );

  my $perform_bacteria_count         = getValue( $def, "perform_bacteria_count", 0 );

  my $blast_top_reads      = getValue( $def, "blast_top_reads" );
  my $blast_unmapped_reads = getValue( $def, "blast_unmapped_reads" );

  my $perform_nonhost_overlap_vis = getValue( $def, "perform_nonhost_overlap_vis", 1 ) && $notMicroRNAOnly;

  my $perform_nonhost_genome_count = getValue( $def, "perform_nonhost_genome_count", 1 ) && $notMicroRNAOnly;

  my $top_read_number = getValue( $def, "top_read_number" );

  my $real_genome_bowtie1_index = $def->{"real_genome_bowtie1_index"};
  my $isHomologyAnalysis        = defined $real_genome_bowtie1_index;

  my $host_genome_dir;
  my $host_genome_suffix = getValue( $def, "host_genome_suffix", "" );
  if ($search_host_genome) {
    $host_genome_dir = create_directory_or_die( $def->{target_dir} . "/host_genome" . $host_genome_suffix );
  }

  my $nonhost_library_dir;
  if ($search_nonhost_library) {
    $nonhost_library_dir = create_directory_or_die( $def->{target_dir} . "/nonhost_library" );
  }

  my $nonhost_genome_dir;
  my @nonhost_genome_groups      = ();
  my @nonhost_genome_group_reads = ();
  my @nonhost_genome_group_names = ();

  my $search_nonhost_genome_custom_group_only = getValue( $def, "search_nonhost_genome_custom_group_only", 0 );
  if ($search_nonhost_genome) {
    $nonhost_genome_dir = create_directory_or_die( $def->{target_dir} . "/nonhost_genome" );
    if ( !$search_nonhost_genome_custom_group_only ) {
      @nonhost_genome_groups      = qw( bacteria_group1 bacteria_group2 fungus_group4 algae_group5 virus_group6 );
      @nonhost_genome_group_reads = qw( bacteria_group1_reads bacteria_group2_reads fungus_group4_reads algae_group5_reads virus_group6_reads );
      @nonhost_genome_group_names = ( "Microbiome Bacteria", "Environment Bacteria", "Fungus", "Algae", "Virus" );
    }

    if ( getValue( $def, "search_nonhost_genome_custom_group", 0 ) ) {
      push( @nonhost_genome_groups,      "custom_group" );
      push( @nonhost_genome_group_reads, "custom_group_reads" );
      push( @nonhost_genome_group_names, getValue( $def, "nonhost_genome_custom_group_name", "Custom" ) );
    }
  }

  my $nonhost_blast_dir;
  if ( ( $blast_unmapped_reads || $search_nonhost_database || $perform_annotate_unmapped_reads ) && ( !$search_nonhost_genome_custom_group_only ) ) {
    $nonhost_blast_dir = create_directory_or_die( $def->{target_dir} . "/final_unmapped" );
  }

  my ( $batchGroups, $batchLayout ) = getBatchGroups($def);
  my $batch_dir;
  if ( defined $batchGroups ) {
    $batch_dir = create_directory_or_die( $def->{target_dir} . "/batch_effects" );
  }

  my $data_visualization_dir = create_directory_or_die( $def->{target_dir} . "/data_visualization" );

  my $perform_tDRmapper = getValue( $def, "perform_tDRmapper", 0 ) && $notMicroRNAOnly;
  if ($perform_tDRmapper) {
    getValue( $def, "tDRmapper" );
    getValue( $def, "tDRmapper_fasta" );
  }

  my $perform_host_tRH_analysis = getValue( $def, "perform_host_tRH_analysis", 0 ) && $notMicroRNAOnly;
  my $perform_host_tRnaFragmentHalves_analysis = getValue( $def, "perform_host_tRnaFragmentHalves_analysis", 0 ) && $notMicroRNAOnly;

  my $R_font_size = 'textSize=9;groupTextSize=' . $def->{table_vis_group_text_size} . ';';

  my $notReadCountPattern = "^.*\/(?!.*?read).*\.count\$";
  my @table_for_correlation = ();
  if ($perform_class_independent_analysis) {
    push @table_for_correlation, ( "identical_sequence_count_table", $notReadCountPattern );
  }

  my @table_for_permanova = ();
  my @table_for_countSum     = ();
  my @table_for_pieSummary   = ();
  my @name_for_pieSummary    = ();
  my @table_for_readSummary  = ();
  my @name_for_readSummary   = ();
  my @table_for_shortReadSource  = ();
  my @name_for_shortReadSource   = ();

  my @name_for_mapPercentage = ( "identical", "dupcount\$" );

  my @reads_for_annoate_unmapped  = ( "identical", "dupcount\$", "cutadapt", ".fastq.short.gz\$" );
  my @files_for_annotate_unmapped = ();
  my @names_for_annotate_unmapped = ();

  if ( $def->{use_least_groups} ) {
    $def->{correlation_rcode} = $def->{correlation_rcode} . "useLeastGroups<-TRUE;";
  }
  else {
    $def->{correlation_rcode} = $def->{correlation_rcode} . "useLeastGroups<-FALSE;";
  }

  if ( $def->{show_label_PCA} ) {
    $def->{correlation_rcode} = $def->{correlation_rcode} . "showLabelInPCA<-TRUE;";
  }
  else {
    $def->{correlation_rcode} = $def->{correlation_rcode} . "showLabelInPCA<-FALSE;";
  }

  my $normalize_by = getValue($def, "normalize_by");

  if ( $def->{correlation_rcode} !~ /totalCountKey/ ) {
    my $correlation_totalCountKey = $normalize_by eq 'None' ? 'None':$normalize_by eq 'TotalReads'?'Reads for Mapping':$normalize_by;
    $def->{correlation_rcode} = $def->{correlation_rcode} . "totalCountKey='$correlation_totalCountKey';";
  }

  if ( $def->{correlation_rcode} !~ /minMedian/ ) {        #set filter parameters
    $def->{correlation_rcode} = $def->{correlation_rcode} . "minMedian=0;minMedianInGroup=1;";
  }

  #print Dumper($config);
  my $groups = $def->{groups};
  if ( !defined $def->{groups_vis_layout} && defined $groups && defined $groups->{".order"} && defined $groups->{".col"} && defined $groups->{".row"} ) {
    $def->{groups_vis_layout} = {
      "Col_Group" => $groups->{".col"},
      "Row_Group" => $groups->{".row"},
      "Groups"    => $groups->{".order"}
    };
  }

  if ( !defined $def->{groups_smallRNA_vis_layout} ) {
    $def->{groups_smallRNA_vis_layout} = $def->{groups_vis_layout};
  }

  my $libraryKey = getValue( $def, "DE_library_key", "TotalReads" );
  my $libraryFile = undef;
  if ( $libraryKey ne "" ) {
    if ($search_host_genome) {
      $libraryFile = [ "bowtie1_genome_1mm_NTA_smallRNA_category", ".Category.Table.csv" ];
    }
    elsif ( defined $config->{fastqc_post_trim_summary} ) {
      $libraryFile = [ "fastqc_post_trim_summary", ".FastQC.reads.tsv\$" ];
      $libraryKey = "Reads";
    }
    elsif ( defined $config->{fastqc_raw_summary} ) {
      $libraryFile = [ "fastqc_raw_summary", ".FastQC.reads.tsv\$" ];
      $libraryKey = "Reads";
    }
    else {
      die( "I don't know where to get library size for key " . $libraryKey );
    }
  }
  else {
    $libraryKey = undef;
  }

  my $hostLibraryStr    = "";
  my $nonhostLibraryStr = "";
  if ( defined $libraryKey ) {
    if ( $libraryKey ne "None" ){
      $hostLibraryStr = $libraryKey;
    }
    if ( $libraryKey eq "TotalReads" ) {
      $nonhostLibraryStr = $hostLibraryStr;
    }
  }

  #print(Dumper($def->{pairs_config}));
  my $do_comparison = (defined $def->{pairs}) || (defined $def->{pairs_config});
  print("do_comparison=" . $do_comparison . "\n");
  if ($do_comparison) {
    my $pairs = $def->{pairs};
    if(not defined $pairs){
      $pairs = $def->{pairs_config};
    }

    my $sampleComparisons;
    if ( defined $pairs->{".order"} ) {
      $sampleComparisons = $pairs->{".order"};
    }
    else {
      my @tmp = sort keys %$pairs;
      $sampleComparisons = \@tmp;
    }

    my $comparisons;
    if ( defined $pairs->{".col"} ) {
      $comparisons = $pairs->{".col"};
    }
    else {
      $comparisons = $sampleComparisons;
    }

    my $hostSmallRNA       = [];
    my $hostSmallRNAFolder = [];
    if(getValue($def, "use_isomiR_in_vis", 1)){
      push( @$hostSmallRNA,       "isomiR" );
      push( @$hostSmallRNAFolder, "miRNA_isomiR" );
    }else{
      push( @$hostSmallRNA,       "miRNA" );
      push( @$hostSmallRNAFolder, "miRNA" );
    }

    if ($notMicroRNAOnly) {
      push( @$hostSmallRNA,       "tDR-anticodon" );
      push( @$hostSmallRNAFolder, "tRNA" );
      if ( $def->{hasSnRNA} ) {
        push( @$hostSmallRNA,       "snDR" );
        push( @$hostSmallRNAFolder, "snRNA" );
      }
      if ( $def->{hasSnoRNA} ) {
        push( @$hostSmallRNA,       "snoDR" );
        push( @$hostSmallRNAFolder, "snoRNA" );
      }
      if ( $def->{hasYRNA} ) {
        push( @$hostSmallRNA,       "yDR" );
        push( @$hostSmallRNAFolder, "yRNA" );
      }
      push( @$hostSmallRNA,       ( "rDR",  "osRNA" ) );
      push( @$hostSmallRNAFolder, ( "rRNA", "otherSmallRNA" ) );

      if ( $def->{hasERV} ) {
        push( @$hostSmallRNA,       "ERV" );
        push( @$hostSmallRNAFolder, "ERV" );
      }
    }
    my $numberOfHostSmallRNA = scalar(@$hostSmallRNA);

    my $DE_task_suffix = getValue( $def, "DE_task_suffix", "" );

    my $numberOfComparison = scalar(@$sampleComparisons);
    if ( !defined $def->{pairs_top_deseq2_vis_layout} ) {
      $def->{pairs_top_deseq2_vis_layout} = {
        "Col_Group" => $comparisons,
        "Row_Group" => [ ("Top 100") x $numberOfComparison ],
        "Groups"    => string_combination( [ ["top100"], [ $nonhostLibraryStr . $DE_task_suffix ], $sampleComparisons ], '_' ),
      };
    }

    if ( !defined $def->{pairs_host_deseq2_vis_layout} ) {
      $def->{pairs_host_deseq2_vis_layout} = {
        "Col_Group" => [ (@$comparisons) x $numberOfHostSmallRNA ],
        "Row_Group" => string_repeat( $hostSmallRNA, $numberOfComparison ),
        "Groups" => string_combination( [ $hostSmallRNAFolder, [ $hostLibraryStr . $DE_task_suffix ], $sampleComparisons ], '_' ),
      };
    }
    if ( !defined $def->{pairs_host_reads_deseq2_vis_layout} ) {
      my $hostSmallRNAReadsFolder = $hostSmallRNAFolder;
      s/miRNA_isomiR/miRNA/ for @{$hostSmallRNAReadsFolder};

      my $hostSmallRNAReads = $hostSmallRNA;
      s/^isomiR$/miRNA/ for @{$hostSmallRNAReads};

      $def->{pairs_host_reads_deseq2_vis_layout} = {
        "Col_Group" => [ (@$comparisons) x $numberOfHostSmallRNA ],
        "Row_Group" => string_repeat( $hostSmallRNAReads, $numberOfComparison ),
        "Groups" => string_combination( [ $hostSmallRNAReadsFolder, [ $hostLibraryStr . $DE_task_suffix ], $sampleComparisons ], '_' ),
      };
    }

    if ( !defined $def->{pairs_host_miRNA_deseq2_vis_layout} ) {
      $def->{pairs_host_miRNA_deseq2_vis_layout} = {
        "Col_Group" => [ (@$comparisons) x 3 ],
        "Row_Group" => string_repeat( [ "isomiR", "NTA", "isomiR NTA" ], $numberOfComparison ),
        "Groups" => string_combination( [ ["miRNA"], [ "isomiR", "NTA", "isomiR_NTA" ], [ $hostLibraryStr . $DE_task_suffix ], $sampleComparisons ], '_' ),
      };
    }

    if ( !defined $def->{pairs_nonHostGroups_deseq2_vis_layout} ) {
      $def->{pairs_nonHostGroups_deseq2_vis_layout} = {
        "Col_Group" => [ (@$comparisons) x scalar(@nonhost_genome_groups) ],
        "Row_Group" => string_repeat( \@nonhost_genome_group_names, $numberOfComparison ),
        "Groups" => string_combination( [ \@nonhost_genome_groups, [ $nonhostLibraryStr . $DE_task_suffix ], $sampleComparisons ], '_' ),
      };
    }

    if ( !defined $def->{pairs_nonHostLibrary_deseq2_vis_layout} ) {
      $def->{pairs_nonHostLibrary_deseq2_vis_layout} = {
        "Col_Group" => [ (@$comparisons) x 5 ],
        "Row_Group" => string_repeat( [ "tDR", "tDR Species", "tDR Amino Acid", "tDR Anticodon", "tDR Reads" ], $numberOfComparison ),
        "Groups" => string_combination( [ ["nonhost_tRNA"], [ "", "species", "type", "anticodon", "reads" ], [ $nonhostLibraryStr . $DE_task_suffix ], $sampleComparisons ], '_' )
      };
    }

    $def->{pure_pairs} = get_pure_pairs( $pairs );
  }

  my $perform_nonhost_tRNA_bacteria_reads = getValue( $def, "perform_nonhost_tRNA_bacteria_reads", 0 );

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

  my $identical_count_ref = [$identical_ref->[0], ".dupcount\$" ];

  my $host_intermediate_dir = getIntermidiateDir($host_genome_dir, $def);

  my $short_reads_table;
  
  if ($search_host_genome) {
    getValue( $def, "coordinate" );

    #1 mismatch search, NTA
    my $hostBowtieTask = "bowtie1_genome_1mm_NTA";
    addBowtie( $config, $def, $individual_ref, $hostBowtieTask, $host_intermediate_dir, $def->{bowtie1_index}, $host_identical_ref, $def->{bowtie1_option_1mm} );

    my $bamSource = $hostBowtieTask;

    if(getValue($def, "is_mix_genome", 0)){
      my $filterMixBam_folder = $host_intermediate_dir . "/" . getNextFolderIndex($def) . "filterMixBam";
      my $host_prefix = getValue($def, "host_prefix");
      $config->{"filterMixBam"} = {
        class                 => "CQS::ProgramWrapperOneToOne",
        perform               => 1,
        target_dir            => $filterMixBam_folder,
        init_command          => "echo sort bam by query name
samtools sort -o __NAME__.name.bam -n __FILE__

echo filter bam by host_prefix $host_prefix ",
        post_command          => "echo sort bam by coordinate
samtools sort -o __NAME__.fixed.bam __NAME__.filtered.bam
rm -f __NAME__.name.bam __NAME__.filtered.bam
mv __NAME__.filtered.txt __NAME__.fixed.txt
",
        option                => "-i __NAME__.name.bam -o __NAME__.filtered.bam --host_prefix " . $host_prefix,
        interpretor           => "python3",
        check_program         => 1,
        program               => "../Alignment/filterMixBam.py",
        source_ref            => [ $bamSource, ".bam\$" ],
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

      $bamSource = ["filterMixBam", ".fixed.bam"];
      push @$individual_ref, ("filterMixBam");

      # my $filterMixBam_summary = "filterMixBam_summary";
      # $config->{$filterMixBam_summary} = {
      #   class         => "CQS::UniqueR",
      #   perform       => 1,
      #   target_dir    => $host_genome_dir . "/" . getNextFolderIndex($def) . "$filterMixBam_summary",
      #   program       => "../Count/count_table.py",
      #   option        => "-p ENS --noheader",
      #   source_arg    => "-i",
      #   source_ref    => $dexseq_count,
      #   output_arg    => "-o",
      #   output_file_ext => ".count",
      #   sh_direct     => 1,
      #   pbs           => {
      #     "nodes"    => "1:ppn=1",
      #     "walltime" => "10",
      #     "mem"      => "20gb"
      #   },
      # };
      # push @$summary, "$dexseq_count_table";
    }

    if ($isHomologyAnalysis) {
      my $realBowtieTask = "bowtie1_real_genome_1mm_NTA";
      addBowtie( $config, $def, $individual_ref, $realBowtieTask, $host_intermediate_dir, $real_genome_bowtie1_index, $host_identical_ref, $def->{bowtie1_option_1mm} );

      my $homologyTask = "bowtie1_genome_1mm_NTA_homology";
      $config->{$homologyTask} = {
        class              => "SmallRNA::FilterIndividualHomologyBAM",
        perform            => 1,
        target_dir         => $host_intermediate_dir . "/$homologyTask",
        option             => "",
        samonly            => 0,
        source_ref         => $hostBowtieTask,
        reference_bams_ref => $realBowtieTask,
        sh_direct          => 1,
        cluster            => $cluster,
        pbs                => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "12",
          "mem"       => "40gb"
        },
      };

      $bamSource = $homologyTask;
      push @$individual_ref, ("$homologyTask");
    }

    my $host_genome = {

      bowtie1_genome_1mm_NTA_smallRNA_count => {
        class           => "CQS::SmallRNACount",
        perform         => 1,
        target_dir      => $host_intermediate_dir . "/bowtie1_genome_1mm_NTA_smallRNA_count",
        option          => $def->{host_smallrnacount_option},
        source_ref      => $bamSource,
        fastq_files_ref => $host_identical_ref,
        seqcount_ref    => $identical_count_ref,
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
    };
    if ( $notMicroRNAOnly and defined $config->{identical_check_cca} ) {
      $host_genome->{bowtie1_genome_1mm_NTA_smallRNA_count}{cca_files_ref} = ["identical_check_cca"];
    }

    push @$individual_ref, ("bowtie1_genome_1mm_NTA_smallRNA_count");

    my $countTask = "bowtie1_genome_1mm_NTA_smallRNA_count";

    my $category_ext = ".Category.Table.csv";
    my $category_other_ext = ".Category1.Barplot.png;.Category2.Barplot.png";
    if (defined $def->{groups}) {
      $category_other_ext = $category_other_ext . ";.Category1.Group.Piechart.png;.Category2.Group.Piechart.png";
    }

    $host_genome = merge_hash_right_precedent(
      $host_genome,
      {
        bowtie1_genome_1mm_NTA_smallRNA_table => {
          class      => "CQS::SmallRNATable",
          perform    => 1,
          target_dir => $host_genome_dir . "/bowtie1_genome_1mm_NTA_smallRNA_table",
          option     => $def->{host_smallrnacounttable_option},
          source_ref => [ $countTask, ".mapped.xml" ],
          prefix     => "smallRNA_1mm_",
          sh_direct  => 1,
          cluster    => $cluster,
          pbs        => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "10",
            "mem"       => "40gb"
          },
        },
        bowtie1_genome_1mm_NTA_smallRNA_info => {
          class      => "CQS::CQSDatatable",
          perform    => 1,
          target_dir => $host_genome_dir . "/bowtie1_genome_1mm_NTA_smallRNA_table",
          option     => "--noheader",
          source_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".info" ],
          prefix     => "smallRNA_1mm_",
          suffix     => ".mapped",
          sh_direct  => 1,
          cluster    => $cluster,
          pbs        => {
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
          output_file_ext           => $category_ext,
          output_other_ext          => $category_other_ext,
          parameterSampleFile1_ref  => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".info" ],
          parameterSampleFile2      => $groups,
          parameterSampleFile2Order => $def->{groups_order},
          parameterSampleFile3      => $def->{groups_smallRNA_vis_layout},
          rCode                     => $R_font_size,
          sh_direct                 => 1,
          pbs                       => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        },
      }
    );
    $config = merge_hash_right_precedent( $config, $host_genome );
    push @$summary_ref, ( "bowtie1_genome_1mm_NTA_smallRNA_table", "bowtie1_genome_1mm_NTA_smallRNA_info", "bowtie1_genome_1mm_NTA_smallRNA_category" );

    my $tRnaAnalysis = {};
    if ($notMicroRNAOnly) {
      $config->{host_genome_tRNA_category} = {
        class                     => "CQS::UniqueR",
        perform                   => 1,
        target_dir                => $data_visualization_dir . "/host_genome_tRNA_category",
        rtemplate                 => "countTableVisFunctions.R,hostTrnaMappingVis.R",
        output_file               => ".tRNAMapping.Result",
        output_file_ext           => ".tRNAType1.Barplot.png;.tRNAType2.Barplot.png",
        parameterSampleFile1Order => $def->{groups_order},
        parameterSampleFile1      => $groups,
        parameterSampleFile2      => $def->{groups_vis_layout},
        parameterFile1_ref        => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.count\$" ],
        parameterFile3_ref        => [ "fastqc_count_vis", ".Reads.csv\$" ],
        rCode                     => 'maxCategory=3;' . $R_font_size,
        sh_direct                 => 1,
        pbs                       => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push @$summary_ref, ("host_genome_tRNA_category");

      if ($perform_host_tRnaFragmentHalves_analysis) {
        $tRnaAnalysis->{"tFragment"} = {
          minLength => 0,
          maxLength => 29
        };
        $tRnaAnalysis->{"tHalves"} = {
          minLength => 30,
          maxLength => 100
        };
      }

      if($perform_host_tRH_analysis){
        $tRnaAnalysis->{"tRH"} = {
          minLength => 30,
          maxLength => 40
        };
      }

      for my $tKey (sort keys %$tRnaAnalysis){
        my $minLength = $tRnaAnalysis->{$tKey}{minLength};
        my $maxLength = $tRnaAnalysis->{$tKey}{maxLength};

        my $tTask = "bowtie1_genome_1mm_NTA_smallRNA_count_" . $tKey;
        my $tTableTask = $tTask . "_table";
        my $tCategoryTask = $tTableTask . "_category";

        $tRnaAnalysis->{$tKey}{"count_section"} = $tTask;
        $tRnaAnalysis->{$tKey}{"table_section"} = $tTableTask;
        $tRnaAnalysis->{$tKey}{"category_section"} = $tCategoryTask;

        $config->{$tTask} = {
          class                 => "CQS::ProgramIndividualWrapper",
          perform               => 1,
          target_dir            => $host_intermediate_dir . "/$tTask",
          option                => "--minLength " . $minLength . " --maxLength " . $maxLength,
          interpretor           => "python3",
          program               => "../SmallRNA/filterTrnaXml.py",
          source_arg            => "-i",
          source_ref            => [ $countTask, ".count.mapped.xml" ],
          output_to_same_folder => 1,
          output_arg            => "-o",
          output_ext            => "." . $tKey . ".count.mapped.xml",
          sh_direct             => 1,
          pbs                   => {
            "email"     => $def->{email},
            "emailType" => $def->{emailType},
            "nodes"     => "1:ppn=1",
            "walltime"  => "10",
            "mem"       => "10gb"
          },
        };
        push @$individual_ref, $tTask;

        $config->{$tTableTask} = {
          class      => "CQS::SmallRNATable",
          perform    => 1,
          target_dir => $host_genome_dir . "/$tTableTask",
          option     => $def->{host_smallrnacounttable_option},
          source_ref => [ $tTask, "." . $tKey . ".count.mapped.xml" ],
          prefix     => "smallRNA_1mm_",
          sh_direct  => 1,
          is_tRH     => 1,
          pbs        => {
            "email"     => $def->{email},
            "emailType" => $def->{emailType},
            "nodes"     => "1:ppn=1",
            "walltime"  => "10",
            "mem"       => "40gb"
          },
        };
        push @$summary_ref, $tTableTask;

        $config->{$tCategoryTask} = {
          class                     => "CQS::UniqueR",
          perform                   => 1,
          target_dir                => $host_genome_dir . "/$tCategoryTask",
          rtemplate                 => "countTableVisFunctions.R,hostTrnaMappingVis.R",
          output_file               => "." . $tKey. ".Mapping.Result",
          output_file_ext           => ".tRNAType1.Barplot.png;.tRNAType2.Barplot.png",
          parameterSampleFile1Order => $def->{groups_order},
          parameterSampleFile1      => $groups,
          parameterSampleFile2      => $def->{groups_vis_layout},
          parameterFile1_ref        => [ $tTableTask, ".tRNA.count\$" ],
          parameterFile3_ref        => [ "fastqc_count_vis", ".Reads.csv\$" ],
          rCode                     => 'maxCategory=3;' . $R_font_size,
          sh_direct                 => 1,
          pbs                       => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        };
        push @$summary_ref, $tCategoryTask;
      }

      # if ($perform_host_tRH_analysis) {
      #   $config->{$tRHTask} = {
      #     class                 => "CQS::ProgramIndividualWrapper",
      #     perform               => 1,
      #     target_dir            => $tRH_folder . "/$tRHTask",
      #     option                => "--minLength 30 --maxLength 40",
      #     interpretor           => "python3",
      #     program               => "../SmallRNA/filterTrnaXml.py",
      #     source_arg            => "-i",
      #     source_ref            => [ $countTask, ".count.mapped.xml" ],
      #     output_to_same_folder => 1,
      #     output_arg            => "-o",
      #     output_ext            => ".tRH.count.mapped.xml",
      #     sh_direct             => 1,
      #     pbs                   => {
      #       "email"     => $def->{email},
      #       "emailType" => $def->{emailType},
      #       "nodes"     => "1:ppn=1",
      #       "walltime"  => "10",
      #       "mem"       => "10gb"
      #     },
      #   };
      #   push @$individual_ref, $tRHTask;

      #   $config->{$tRHTableTask} = {
      #     class      => "CQS::SmallRNATable",
      #     perform    => 1,
      #     target_dir => $tRH_folder . "/$tRHTableTask",
      #     option     => $def->{host_smallrnacounttable_option},
      #     source_ref => [ $tRHTask, ".tRH.count.mapped.xml" ],
      #     prefix     => "smallRNA_1mm_",
      #     sh_direct  => 1,
      #     is_tRH     => 1,
      #     cluster    => $cluster,
      #     pbs        => {
      #       "email"     => $def->{email},
      #       "emailType" => $def->{emailType},
      #       "nodes"     => "1:ppn=1",
      #       "walltime"  => "10",
      #       "mem"       => "40gb"
      #     },
      #   };

      #   push @$summary_ref, $tRHTableTask;

      #   $config->{$tRHCategory} = {
      #     class                     => "CQS::UniqueR",
      #     perform                   => 1,
      #     target_dir                => $tRH_folder . "/host_genome_tRH_category",
      #     rtemplate                 => "countTableVisFunctions.R,hostTrnaMappingVis.R",
      #     output_file               => ".tRHMapping.Result",
      #     output_file_ext           => ".tRNAType2.Barplot.png",
      #     parameterSampleFile1Order => $def->{groups_order},
      #     parameterSampleFile1      => $groups,
      #     parameterSampleFile2      => $def->{groups_vis_layout},
      #     parameterFile1_ref        => [ $tRHTableTask, ".tRNA.count\$" ],
      #     parameterFile3_ref        => [ "fastqc_count_vis", ".Reads.csv\$" ],
      #     rCode                     => 'maxCategory=3;' . $R_font_size,
      #     sh_direct                 => 1,
      #     pbs                       => {
      #       "email"     => $def->{email},
      #       "emailType" => $def->{emailType},
      #       "nodes"     => "1:ppn=1",
      #       "walltime"  => "1",
      #       "mem"       => "10gb"
      #     },
      #   };
      #   push @$summary_ref, $tRHCategory;
      # }

      if ( $def->{perform_host_tRNA_start_position} && $def->{groups}) {
        my $tTask = "host_genome_tRNA_start_position_vis";
        if ( !defined $def->{tRNA_vis_group} ) {
          $def->{tRNA_vis_group} = $groups;
        }
        addPositionVis(
          $config, $def,
          $summary_ref,
          $tTask,
          $data_visualization_dir,
          {
            target_dir         => $data_visualization_dir . "/" . $tTask,
            output_file        => ".tRNAStartPositionVis",
            output_file_ext    => ".barplot.png",
            rtemplate          => "tRnaStartPositionVis.R",
            parameterFile1_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.count.startpos\$" ],
          }
        );

        for my $tKey (sort keys %$tRnaAnalysis){
          my $tTableTask = $tRnaAnalysis->{$tKey}{"table_section"};

          my $tStartPositionTask = $tTableTask . "_start_position_vis";
          addPositionVis(
            $config, $def,
            $summary_ref,
            $tStartPositionTask,
            $data_visualization_dir,
            {
              target_dir         => $data_visualization_dir . "/" . $tStartPositionTask,
              output_file        => ".${tKey}.StartPositionVis",
              output_file_ext    => ".barplot.png",
              rtemplate          => "tRnaStartPositionVis.R",
              parameterFile1_ref => [ $tTableTask, ".tRNA.count.startpos\$" ],
            }
          );
        }
      }

      if ( getValue( $def, "perform_host_rRNA_coverage" ) ) {
        my $visualizationTask = "host_genome_rRNA_position_vis";
        my $folder            = $data_visualization_dir . "/" . $visualizationTask;
        $config->{$visualizationTask} = {
          class                    => "CQS::ProgramWrapper",
          perform                  => 1,
          target_dir               => $folder,
          interpretor              => "python3",
          program                  => "../SmallRNA/rRNAHostCoverage.py",
          parameterSampleFile1_arg => "-i",
          parameterSampleFile1_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".mapped.xml" ],
          output_arg               => "-o",
          output_ext               => ".position",
          sh_direct                => 1,
          pbs                      => {
            "email"     => $def->{email},
            "emailType" => $def->{emailType},
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "10gb"
          },
        };
        push( @$summary_ref, $visualizationTask );
      }
    }

    if ( defined $def->{host_xml2bam} && $def->{host_xml2bam} ) {
      $config->{bowtie1_genome_xml2bam} = {
        class         => "SmallRNA::HostXmlToBam",
        perform       => 1,
        target_dir    => $host_intermediate_dir . "/bowtie1_genome_xml2bam",
        source_ref    => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".mapped.xml" ],
        bam_files_ref => [ $bamSource, ".bam" ],
        sh_direct     => 1,
        pbs           => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$individual_ref, "bowtie1_genome_xml2bam" );

      if ( getValue( $def, "host_bamplot" ) ) {
        my $plot_gff = getValue( $def, "host_bamplot_gff" );

        # "-g HG19 -y uniform -r"
        my $bamplot_option = getValue( $def, "host_bamplot_option" );
        my $plotgroups = $def->{plotgroups};
        if ( !defined $plotgroups ) {
          my $files         = $def->{files};
          my @sortedSamples = sort keys %$files;
          $plotgroups = { $def->{task_name} => \@sortedSamples };
        }

        $config->{"plotgroups"}   = $plotgroups;
        $config->{"host_bamplot"} = {
          class              => "Visualization::Bamplot",
          perform            => 1,
          target_dir         => "${host_genome_dir}/host_bamplot",
          option             => $bamplot_option,
          source_ref         => "bowtie1_genome_xml2bam",
          groups_ref         => "plotgroups",
          gff_file           => $plot_gff,
          is_rainbow_color   => 0,
          is_draw_individual => 0,
          is_single_pdf      => 1,
          sh_direct          => 1,
          pbs                => {
            "email"    => $def->{email},
            "nodes"    => "1:ppn=1",
            "walltime" => "1",
            "mem"      => "10gb"
          },
        };
        push @$summary_ref, ("host_bamplot");
      }
    }

    push( @name_for_mapPercentage,      "bowtie1_genome_1mm_NTA_smallRNA_count", "count.mapped.xml" );
    push( @files_for_annotate_unmapped, "bowtie1_genome_1mm_NTA_smallRNA_count", "count.mapped.xml" );
    push( @names_for_annotate_unmapped, "smallRNA" );

    if ( $notMicroRNAOnly && $def->{has_NTA} && $def->{consider_tRNA_NTA} ) {
      $host_genome->{"bowtie1_genome_1mm_NTA_smallRNA_count"}{"cca_file_ref"} = "identical_check_cca";
    }

    push @table_for_pieSummary, ( "bowtie1_genome_1mm_NTA_smallRNA_count", ".count\$" );
    push @name_for_pieSummary, "Host Small RNA";

    push @table_for_shortReadSource, ( "bowtie1_genome_1mm_NTA_smallRNA_table", $task_name . ".read.count\$");
    push @name_for_shortReadSource, ( "host smallRNA");

    push @name_for_readSummary, (
      "Host miRNA",    #miRNA
    );
    push @table_for_readSummary, (
      "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.read.count\$",    #miRNA
    );
    push @table_for_countSum, (
      "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.count\$",         #miRNA
    );
    push @table_for_correlation, (
      "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.count\$",           #miRNA
      "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.isomiR.count\$",    #miRNA isomiR
    );
    if ( $def->{read_correlation} ) {
      push @table_for_correlation, (
        "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.read.count\$",    #miRNA
      );
    }

    if ($notMicroRNAOnly) {
      push @name_for_readSummary, (
        "Host tRNA"                                                        #tRNA
      );
      push @table_for_readSummary, (
        "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.read.count\$"      #tRNA
      );
      push @table_for_countSum, (
        "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.count\$"           #tRNA
      );
      push @table_for_correlation, (
        "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.count\$"           #tRNA
      );

      for my $tKey (sort keys %$tRnaAnalysis){
        my $tTableTask = $tRnaAnalysis->{$tKey}{table_section};

        push @table_for_correlation, (
          $tTableTask, ".tRNA.count\$"                                   #tRNA
        );
      }

      if ( $def->{read_correlation} ) {
        push @table_for_correlation, (
          "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.read.count\$"    #tRNA
        );
        for my $tKey (sort keys %$tRnaAnalysis){
          my $tTableTask = $tRnaAnalysis->{$tKey}{table_section};

          push @table_for_correlation, (
            $tTableTask, ".tRNA.read.count\$"                                   #tRNA
          );
        }
      }
      if ( $def->{hasYRNA} ) {
        push @name_for_readSummary, "Host yRNA";
        push @table_for_countSum,    ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".yRNA.count\$" );
        push @table_for_readSummary, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".yRNA.read.count\$" );
        push @table_for_correlation, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".yRNA.count\$" );
        if ( $def->{read_correlation} ) {
          push @table_for_correlation, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".yRNA.read.count\$", );
        }
      }
      if ( $def->{hasSnRNA} ) {
        push @name_for_readSummary, "Host snRNA";
        push @table_for_countSum,    ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".snRNA.count\$" );
        push @table_for_readSummary, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".snRNA.read.count\$" );
        push @table_for_correlation, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".snRNA.count\$" );
        if ( $def->{read_correlation} ) {
          push @table_for_correlation, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".snRNA.read.count\$", );
        }
      }
      if ( $def->{hasSnoRNA} ) {
        push @name_for_readSummary, "Host snoRNA";
        push @table_for_countSum,    ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".snoRNA.count\$" );
        push @table_for_readSummary, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".snoRNA.read.count\$" );
        push @table_for_correlation, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".snoRNA.count\$" );
        if ( $def->{read_correlation} ) {
          push @table_for_correlation, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".snoRNA.read.count\$", );
        }
      }
      push @name_for_readSummary, (
        "Host rRNA",               #rRNA
        "Host other small RNA",    #other
      );
      push @table_for_readSummary, (
        "bowtie1_genome_1mm_NTA_smallRNA_table", ".rRNA.read.count\$",    #rRNA
        "bowtie1_genome_1mm_NTA_smallRNA_table", ".other.read.count\$"    #other
      );
      push @table_for_countSum, (
        "bowtie1_genome_1mm_NTA_smallRNA_table", ".rRNA.count\$",         #rRNA
        "bowtie1_genome_1mm_NTA_smallRNA_table", ".other.count\$"         #other
      );
      push @table_for_correlation, (
        "bowtie1_genome_1mm_NTA_smallRNA_table", ".rRNA.count\$",         #rRNA
        "bowtie1_genome_1mm_NTA_smallRNA_table", ".other.count\$"         #other
      );

      if ( $def->{hasERV} ) {
        push @name_for_readSummary, "Host ERV";
        push @table_for_countSum,    ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".ERV.count\$" );
        push @table_for_readSummary, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".ERV.read.count\$" );
        push @table_for_correlation, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".ERV.count\$" );
        if ( $def->{read_correlation} ) {
          push @table_for_correlation, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".ERV.read.count\$", );
        }
      }

      if ( $def->{read_correlation} ) {
        push @table_for_correlation, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".rRNA.read.count\$", "bowtie1_genome_1mm_NTA_smallRNA_table", ".other.read.count\$", );
      }
    }

    ##
    if ($do_comparison) {
      my @visual_source       = ();
      my @visual_source_reads = ();

      #miRNA
      addDEseq2( $config, $def, $summary_ref, "miRNA", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.count\$" ], $host_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
      addDEseq2( $config, $def, $summary_ref, "miRNA_NTA", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.NTA.count\$" ],
        $host_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
      addDEseq2( $config, $def, $summary_ref, "miRNA_NTA_base", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.NTA.base.count\$" ],
        $host_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );

      addDEseq2( $config, $def, $summary_ref, "miRNA_isomiR", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.isomiR.count\$" ],
        $host_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );

      if(getValue($def, "use_isomiR_in_vis", 1)){
        push @visual_source, "miRNA_isomiR";
      }else{
        push @visual_source, "miRNA";
      }

      addDEseq2( $config, $def, $summary_ref, "miRNA_isomiR_NTA", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.isomiR_NTA.count\$" ],
        $host_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
      addDEseq2( $config, $def, $summary_ref, "miRNA_reads", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.read.count\$" ],
        $host_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
      addDeseq2Visualization( $config, $def, $summary_ref, "host_genome_miRNA", [ "miRNA_isomiR", "miRNA_NTA", "miRNA_isomiR_NTA" ],
        $data_visualization_dir, "pairs_host_miRNA_deseq2_vis_layout", $libraryKey );
      push @visual_source_reads, "miRNA_reads";

      if ($notMicroRNAOnly) {
        my @biotypes = ();

        addDEseq2( $config, $def, $summary_ref, "tRNA_aminoacid", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.aminoacid.count\$" ],
          $host_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );

        push(@biotypes, "tRNA");

        if ( $def->{hasYRNA} ) {
          push(@biotypes, "yRNA");
        }

        if ( $def->{hasSnRNA} ) {
          push(@biotypes, "snRNA");
        }

        if ( $def->{hasSnoRNA} ) {
          push(@biotypes, "snoRNA");
        }

        push(@biotypes, "rRNA");

        push(@biotypes, "other");

        if ( $def->{hasERV} ) {
          push(@biotypes, "ERV");
        }

        #print(@biotypes);

        for my $biotype (@biotypes){
          my $biotype_name = $biotype eq "other" ? "otherSmallRNA" : $biotype;

          $deseq2Task = addDEseq2( $config, $def, $summary_ref, $biotype_name, [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".${biotype}.count\$" ],
            $host_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
          push( @visual_source, $biotype_name );

          addDEseq2( $config, $def, $summary_ref, "${biotype_name}_reads", [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".${biotype}.read.count\$" ],
            $host_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
          push @visual_source_reads, "${biotype_name}_reads";
        }
      }

      #host genome smallRNA visualization
      addDeseq2Visualization( $config, $def, $summary_ref, "host_genome",       \@visual_source,       $data_visualization_dir, "pairs_host_deseq2_vis_layout",       $libraryKey );
      addDeseq2Visualization( $config, $def, $summary_ref, "host_genome_reads", \@visual_source_reads, $data_visualization_dir, "pairs_host_reads_deseq2_vis_layout", $libraryKey );

      for my $tKey (sort keys %$tRnaAnalysis){
        my $tTableTask = $tRnaAnalysis->{$tKey}{table_section};
        addDEseq2( $config, $def, $summary_ref, $tKey,           [ $tTableTask, ".tRNA.count\$" ],           $host_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
        addDEseq2( $config, $def, $summary_ref, "${tKey}_reads",     [ $tTableTask, ".tRNA.read.count\$" ],      $host_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
        addDEseq2( $config, $def, $summary_ref, "${tKey}_aminoacid", [ $tTableTask, ".tRNA.aminoacid.count\$" ], $host_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
      }
    }

    if ( $do_comparison || defined $groups || defined $def->{tRNA_vis_group} ) {
      addPositionVis(
        $config, $def,
        $summary_ref,
        "host_genome_miRNA_PositionVis",
        $data_visualization_dir,
        {
          output_file        => ".miRNAPositionVis",
          parameterFile1_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.count.position\$" ],
        }
      );

      if ($notMicroRNAOnly) {
        my $trna_sig_result;
        if ( !defined $def->{tRNA_vis_group} ) {
          $def->{tRNA_vis_group} = $groups;
        }

        if ($do_comparison) {
          $trna_sig_result = [ "deseq2_tRNA", "_DESeq2_sig.csv\$" ];
        }

        addPositionVis(
          $config, $def,
          $summary_ref,
          "host_genome_tRNA_PositionVis",
          $data_visualization_dir,
          {
            output_file        => ".tRNAPositionVis",
            parameterFile1_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.aminoacid.count.position\$" ],

            #        parameterSampleFile3_ref => $trna_sig_result,
          }
        );
        addPositionVis(
          $config, $def,
          $summary_ref,
          "host_genome_tRNA_PositionVis_anticodon",
          $data_visualization_dir,
          {
            output_file        => ".tRNAAnticodonPositionVis",
            parameterFile1_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.count.position\$" ],

            #        parameterSampleFile3_ref => $trna_sig_result,
          }
        );
        if ( $def->{hasYRNA} ) {
          addPositionVis(
            $config, $def,
            $summary_ref,
            "host_genome_yRNA_PositionVis",
            $data_visualization_dir,
            {
              output_file        => ".yRNAPositionVis",
              parameterFile1_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".yRNA.count.position\$" ],
            }
          );
        }
        if ( $def->{hasSnRNA} ) {
          addPositionVis(
            $config, $def,
            $summary_ref,
            "host_genome_snRNA_PositionVis",
            $data_visualization_dir,
            {
              output_file        => ".snRNAPositionVis",
              parameterFile1_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".snRNA.count.position\$" ],
            }
          );
        }
        if ( $def->{hasSnoRNA} ) {
          addPositionVis(
            $config, $def,
            $summary_ref,
            "host_genome_snoRNA_PositionVis",
            $data_visualization_dir,
            {
              output_file        => ".snoRNAPositionVis",
              parameterFile1_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_table", ".snoRNA.count.position\$" ],
            }
          );
        }

        for my $tKey (sort keys %$tRnaAnalysis){
          my $tTableTask = $tRnaAnalysis->{$tKey}{table_section};
          addPositionVis(
            $config, $def,
            $summary_ref,
            "host_genome_${tKey}_PositionVis",
            $data_visualization_dir,
            {
              output_file        => ".${tKey}PositionVis",
              parameterFile1_ref => [ $tTableTask, ".tRNA.aminoacid.count.position\$" ],
            }
          );
          addPositionVis(
            $config, $def,
            $summary_ref,
            "host_genome_${tKey}_PositionVis_anticodon",
            $data_visualization_dir,
            {
              output_file        => ".${tKey}AnticodonPositionVis",
              parameterFile1_ref => [ $tTableTask, ".tRNA.count.position\$" ],
            }
          );
        }
      }
    }

    if ( $search_nonhost_database || $blast_unmapped_reads || $def->{perform_host_length_dist_category} || $def->{perform_host_genome_reads_deseq2}  || $def->{perform_short_reads_deseq2} || $def->{perform_short_reads_source}) {
      my $readClass;
      my $readTask;
      if ( $def->{host_remove_all_mapped_reads} ) {
        $readClass = "Samtools::MappedReadNames";
        $readTask  = "bowtie1_genome_1mm_NTA_read_names";
      }
      else {
        $readClass = "Samtools::PerfectMappedReadNames";
        $readTask  = "bowtie1_genome_1mm_NTA_pmnames";
      }
      my $unmapped_reads = {

        #matched reads with host genome
        $readTask => {
          class      => $readClass,
          perform    => 1,
          target_dir => $host_intermediate_dir . "/" . $readTask,
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
          target_dir  => $host_intermediate_dir . "/bowtie1_genome_unmapped_reads",
          perlFile    => "unmappedReadsToFastq.pl",
          source_ref  => $identical_ref,
          source2_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_count", ".mapped.xml" ],
          source3_ref => [$readTask],
          output_ext  => "_clipped_identical.unmapped.fastq.gz",
          output_other_ext =>
"_clipped_identical.unmapped.fastq.dupcount,_clipped_identical.mappedToHostGenome.fastq.gz,_clipped_identical.mappedToHostGenome.fastq.dupcount,_clipped_identical.short.fastq.gz,_clipped_identical.short.fastq.dupcount,_clipped_identical.unmapped.fastq.gz.info",
          sh_direct => 1,
          pbs       => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => getValue($def, "bowtie1_genome_unmapped_reads_mem", "10gb")
          },
        },
        bowtie1_genome_host_reads_table => {
          class      => "CQS::CQSDatatable",
          perform    => 1,
          target_dir => $host_genome_dir . "/bowtie1_genome_host_reads_table",
          source_ref => [ "bowtie1_genome_unmapped_reads", ".mappedToHostGenome.fastq.dupcount\$" ],
          option     => "-k 2 -v 1 --fillMissingWithZero",
          sh_direct  => 1,
          pbs        => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "1",
            "mem"       => "40gb"
          },
        },
      };
      $config = merge_hash_right_precedent( $config, $unmapped_reads );
      push @table_for_shortReadSource, ( "bowtie1_genome_host_reads_table", ".count\$");
      push @name_for_shortReadSource, ( "host genome");

      push( @name_for_mapPercentage,      "bowtie1_genome_unmapped_reads", ".mappedToHostGenome.fastq.dupcount\$" );
      push( @files_for_annotate_unmapped, "bowtie1_genome_unmapped_reads", ".mappedToHostGenome.fastq.dupcount\$" );
      push( @names_for_annotate_unmapped, "host_genome" );

      push @$individual_ref, ( $readTask, "bowtie1_genome_unmapped_reads" );
      push @$summary_ref, ("bowtie1_genome_host_reads_table");
      push @table_for_pieSummary,
        ( "bowtie1_genome_unmapped_reads", ".mappedToHostGenome.fastq.dupcount", "bowtie1_genome_unmapped_reads", ".short.fastq.dupcount", "bowtie1_genome_unmapped_reads",
        ".unmapped.fastq.dupcount" );
      push @name_for_pieSummary, ( "Mapped to Host Genome", "Too Short for Mapping", "Unmapped In Host" );
      push @table_for_readSummary, ( "bowtie1_genome_host_reads_table", ".count\$" );
      push @name_for_readSummary, ("Host Genome");
      $identical_ref       = [ "bowtie1_genome_unmapped_reads", ".unmapped.fastq.gz\$" ];
      $identical_count_ref = [ "bowtie1_genome_unmapped_reads", ".unmapped.fastq.dupcount\$" ];

      if ( $do_comparison and $def->{perform_host_genome_reads_deseq2} ) {
        addDEseq2( $config, $def, $summary_ref, "bowtie1_genome_host_reads", [ "bowtie1_genome_host_reads_table", ".count\$" ],
          $host_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
      }

      if($def->{perform_short_reads_deseq2} || $def->{perform_short_reads_source}){
        $short_reads_table = "bowtie1_genome_short_reads_table";
        $config->{$short_reads_table} = {
          class                 => "CQS::ProgramWrapper",
          perform               => 1,
          target_dir            => "$host_genome_dir/bowtie1_genome_short_reads_table",
          option                => "",
          interpretor           => "python3",
          check_program         => 1,
          program               => "../SmallRNA/shortReadTable.py",
          source_ref            => [ "bowtie1_genome_unmapped_reads", "_clipped_identical.short.fastq.dupcount" ],
          source_arg            => "-i",
          output_arg            => "-o",
          output_file_prefix    => "",
          output_file_ext       => ".count.txt",
          output_other_ext      => "",
          sh_direct             => 1,
          pbs                   => {
            "nodes"     => "1:ppn=1",
            "walltime"  => "10",
            "mem"       => "10gb"
          },
        };
        push @$summary_ref, ("bowtie1_genome_short_reads_table");

        if ( $def->{perform_short_reads_deseq2} && $do_comparison ) {
          addDEseq2( $config, $def, $summary_ref, "short_reads", [ $short_reads_table, ".count.txt\$" ],
            $host_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
        }
      }
    }

    if ( $def->{perform_host_length_dist_category} ) {
      my @length_dist_count = ();
      my @length_dist_names = ();

      push @length_dist_names, ( "miRNA", "tDR", "rDR" );
      push @length_dist_count, (
        "bowtie1_genome_1mm_NTA_smallRNA_table", ".miRNA.read.count\$",
        "bowtie1_genome_1mm_NTA_smallRNA_table", ".tRNA.read.count\$",
        "bowtie1_genome_1mm_NTA_smallRNA_table", ".rRNA.read.count\$"     #rRNA
      );

      if ( $def->{hasYRNA} ) {
        push @length_dist_names, "yDR";
        push @length_dist_count, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".yRNA.read.count\$" );
      }

      if ( $def->{hasSnRNA} ) {
        push @length_dist_names, "snDR";
        push @length_dist_count, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".snRNA.read.count\$" );
      }

      if ( $def->{hasSnoRNA} ) {
        push @length_dist_names, "snoDR";
        push @length_dist_count, ( "bowtie1_genome_1mm_NTA_smallRNA_table", ".snoRNA.read.count\$" );
      }

      push @length_dist_names, ( "osDR", "Genome", "fastq_len", "category" );
      push @length_dist_count, (
        "bowtie1_genome_1mm_NTA_smallRNA_table",    ".other.read.count\$",         #other
        "bowtie1_genome_host_reads_table",          ".count\$",
        "fastq_len_vis",                            ".lengthDistribution.csv\$",
        "bowtie1_genome_1mm_NTA_smallRNA_category", "Category.Table.csv\$"
      );

      $config->{host_length_dist_category} = {
        class                     => "CQS::UniqueR",
        perform                   => 1,
        target_dir                => $data_visualization_dir . "/host_length_dist_category",
        rtemplate                 => "countTableVisFunctions.R,../SmallRNA/lengthDistributionStackedBarplot.R",
        output_file               => ".length",
        output_file_ext           => ".pdf;.png",
        parameterSampleFile1_ref  => \@length_dist_count,
        parameterSampleFile1Names => \@length_dist_names,
        parameterSampleFile2      => $def->{groups},
        sh_direct                 => 1,
        rCode                     => '' . $R_font_size,
        pbs                       => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push @$summary_ref, ("host_length_dist_category");
    }
  }

  if ( defined $libraryKey && $libraryKey ne "TotalReads" ) {
    $libraryFile = undef;
    $libraryKey  = undef;
  }

  if ( $do_comparison and $perform_class_independent_analysis ) {
    my $taskKey = "top${top_read_number}";

    addDEseq2( $config, $def, $summary_ref, "${taskKey}_reads", [ "identical_sequence_count_table", ".read.count\$" ], $class_independent_dir, $DE_min_median_read_top, $libraryFile, $libraryKey );
    addDeseq2Visualization( $config, $def, $summary_ref, "${taskKey}_reads", ["${taskKey}_reads"], $data_visualization_dir, "pairs_top_deseq2_vis_layout", $libraryKey );

    addDEseq2( $config, $def, $summary_ref, "${taskKey}_contigs", [ "identical_sequence_count_table", ".count\$" ], $class_independent_dir, $DE_min_median_read_top, $libraryFile, $libraryKey );
    addDeseq2Visualization( $config, $def, $summary_ref, "${taskKey}_contigs", ["${taskKey}_contigs"], $data_visualization_dir, "pairs_top_deseq2_vis_layout", $libraryKey );

    addDEseq2( $config, $def, $summary_ref, "${taskKey}_minicontigs", [ "identical_sequence_count_table", ".minicontig.count\$" ],
      $class_independent_dir, $DE_min_median_read_top, $libraryFile, $libraryKey );
    addDeseq2Visualization( $config, $def, $summary_ref, "${taskKey}_minicontigs", ["${taskKey}_minicontigs"], $data_visualization_dir, "pairs_top_deseq2_vis_layout", $libraryKey );
  }

  if ( $search_nonhost_database && getValue( $def, "search_combined_nonhost" ) ) {

    #Mapping host genome reads to non-host databases
    addNonhostDatabase(
      $config, $def, $individual_ref, $summary_ref, "HostGenomeReads_NonHost_pm", $nonhost_library_dir,    #general option
      $def->{bowtie1_all_nonHost_index}, [ "bowtie1_genome_unmapped_reads", ".mappedToHostGenome.fastq.gz" ],    #bowtie option
      $def->{smallrnacount_option} . ' --keepChrInName --categoryMapFile ' . $def->{all_nonHost_map},            #count option
      $def->{nonhost_table_option},                                                                              #table option
      $identical_count_ref
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

  my $nonhost_count = [];
  my $nonhost_genome_read_count = [];
  my $bacteria_read_count = [];

  if (getValue($def, "search_nonhost_genome", 1) && getValue($def, "search_refseq_bacteria", 0)) {
    my $refseq_bacteria_bowtie = "refseq_bacteria_bowtie";
    $config->{"refseq_bacteria_bowtie_index"} = getValue($def, "refseq_bacteria_bowtie_index");
    $config->{$refseq_bacteria_bowtie} = {
      class => "CQS::ProgramWrapperOneToManyScatter",
      target_dir => "$host_intermediate_dir/$refseq_bacteria_bowtie",
      interpretor => "",
      program => "",
      check_program => 0,
      option => "
rm -f __NAME_____SCATTERNAME__.failed

bowtie -a --best -m 10000 --strata -v 0 -p 8 -x __SCATTER__ __FILE__ __OUTPUT__.tmp

status=\$?
if [[ \$status -ne 0 ]]; then
  touch __NAME_____SCATTERNAME__.failed
  rm __OUTPUT__.tmp
else
  cut -f1-4 __OUTPUT__.tmp > __OUTPUT__
  gzip __OUTPUT__
  rm __OUTPUT__.tmp
fi
",
      source_arg => "-x",
      source_ref => $identical_ref,
      scatter_arg => "",
      scatter_ref => "refseq_bacteria_bowtie_index",
      output_arg => "",
      output_file_prefix => ".txt",
      output_file_ext => ".txt.gz",
      output_to_same_folder => 0,
      can_result_be_empty_file => 0,
      use_tmp_folder => getValue($def, "use_tmp_folder", 0),
      sh_direct   => 0,
      pbs => {
        "nodes"     => "1:ppn=" . getValue($def, "${refseq_bacteria_bowtie}_nodes", "8"),,
        "walltime"  => getValue($def, "${refseq_bacteria_bowtie}_walltime", "4"),
        "mem"       => getValue($def, "${refseq_bacteria_bowtie}_mem", "20gb"),
      },
    };

    push( @$individual_ref, $refseq_bacteria_bowtie );

    my $refseq_bacteria_bowtie_count = "refseq_bacteria_bowtie_count";
    $config->{$refseq_bacteria_bowtie_count} = {
      class => "CQS::ProgramWrapperManyToOneGather",
      target_dir => "$host_intermediate_dir/$refseq_bacteria_bowtie_count",
      interpretor => "",
      program => getValue($def, "spcount", "spcount"),
      check_program => 0,
      option => "bowtie_count --species_column " . getValue($def, "species_column", "species"),
      source_arg => "-i",
      source_ref => $identical_ref,
      parameterSampleFile2_arg => "-c",
      parameterSampleFile2_ref => $identical_count_ref,
      sample_scatter_ref => $refseq_bacteria_bowtie,
      scatter_ref => "refseq_bacteria_bowtie_index",
      parameterFile1_arg => "-s",
      parameterFile1 => getValue($def, "refseq_bacteria_species"),
      output_arg => "-o",
      output_file_prefix => ".txt.gz",
      output_file_ext => ".txt.gz",
      output_to_same_folder => 1,
      can_result_be_empty_file => 0,
      #no_docker => 1,
      use_tmp_folder => getValue($def, "use_tmp_folder", 0),
      sh_direct   => 0,
      pbs => {
        "nodes"     => "1:ppn=1",
        "walltime"  => getValue($def, "${refseq_bacteria_bowtie_count}_walltime", "2"),
        "mem"       => getValue($def, "${refseq_bacteria_bowtie_count}_mem", "20gb"),
      },
    };

    push( @$individual_ref, $refseq_bacteria_bowtie_count );

    my $categories = ["species", "genus", "family", "order", "class", "phylum" ];
    my $file_exts = [];
    for my $cat (@$categories){
      push(@$file_exts, ".${cat}.query.count");
      push(@$file_exts, ".${cat}.estimated.count");
      push(@$file_exts, ".${cat}.aggregated.count");
    }
    my $file_ext_str = ".read.count,.tree.count," . join(",", @$file_exts);

    my $refseq_bacteria_table = "refseq_bacteria_table";
    $config->{$refseq_bacteria_table} = {
      class => "CQS::ProgramWrapper",
      target_dir => "$nonhost_genome_dir/$refseq_bacteria_table",
      interpretor => "",
      program => getValue($def, "spcount", "spcount"),
      check_program => 0,
      option => "count_table -o __NAME__ --species_column " . getValue($def, "species_column", "species"),
      source_arg => "-i",
      source_ref => $refseq_bacteria_bowtie_count,
      parameterFile1_arg => "-s",
      parameterFile1 => getValue($def, "refseq_bacteria_species"),
      parameterFile2_arg => "-t",
      parameterFile2 => getValue($def, "refseq_taxonomy"),
      output_arg => "-o",
      output_file_ext => $file_ext_str,
      #no_docker => 1,
      sh_direct   => 1,
      pbs => {
        "nodes"     => "1:ppn=1",
        "walltime"  => getValue($def, "${refseq_bacteria_table}_walltime", "24"),
        "mem"       => getValue($def, "${refseq_bacteria_table}_mem", "100gb"),
      },
    };

    push( @$summary_ref, $refseq_bacteria_table );  
    push @table_for_correlation, ( $refseq_bacteria_table, ".estimated.count\$" );
    push @table_for_correlation, ( $refseq_bacteria_table, ".aggregated.count\$" );

    push @table_for_countSum,    ( $refseq_bacteria_table, ".phylum.estimated.count\$" );

    push @table_for_readSummary, ( $refseq_bacteria_table, '.read.count$' );
    push @name_for_readSummary,  ( "Refseq Bacteria" );

    push @$nonhost_count, ( $refseq_bacteria_table, '.read.count$' );
    push @$nonhost_genome_read_count, ( $refseq_bacteria_table, '.read.count$' );
    push @$bacteria_read_count, ( $refseq_bacteria_table, '.read.count$' );
  
    my @file_exts = ();

    my $task_name = $def->{task_name};
    my $files = $def->{files};
    my $groups = $def->{groups};
    for my $sample (sort keys %$files){
      push(@file_exts, ".${sample}.html");
    }
    for my $gname (sort keys %$groups){
      push(@file_exts, ".${gname}.html");
    }

    my $krona_count_table = getValue($def, "krona_count_table", "estimated");
    my $refseq_bacteria_krona = "refseq_bacteria_krona_" . $krona_count_table;
    my $krona_ref;
    if ($krona_count_table eq "estimated"){
      $krona_ref = [$refseq_bacteria_table, ".species.estimated.count"];
    }elsif($krona_count_table eq "aggregated"){
      $krona_ref = [$refseq_bacteria_table, ".tree.count"];
    }else{
      die "krona_count_table should be either estimated or aggregated, now is " . $krona_count_table;
    }

    $config->{$refseq_bacteria_krona} = {
      class => "CQS::ProgramWrapper",
      target_dir => "$data_visualization_dir/$refseq_bacteria_krona",
      interpretor => "",
      program => getValue($def, "spcount", "spcount"),
      check_program => 0,
      option => "krona -o __NAME__ ",
      post_command => "rm -rf *.html.files .cache .config",
      parameterSampleFile1_arg => "-g",
      parameterSampleFile1 => $groups,
      parameterFile1_arg => "-i",
      parameterFile1_ref => $krona_ref,
      parameterFile2_arg => "-t",
      parameterFile2     => getValue($def, "krona_taxonomy_folder"),
      output_arg => "-o",
      output_file_ext => join(",", @file_exts),
      #no_docker => 1,
      sh_direct   => 1,
      pbs => {
        "nodes"     => "1:ppn=1",
        "walltime"  => getValue($def, "${refseq_bacteria_krona}_walltime", "10"),
        "mem"       => getValue($def, "${refseq_bacteria_krona}_mem", "20gb"),
      },
    };
    push( @$summary_ref, $refseq_bacteria_krona );  

    if ($do_comparison) {
      for my $cat (@$categories){
        # addDEseq2( $config, $def, $summary_ref, "refseq_bacteria_${cat}_query", [ $refseq_bacteria_table, ".${cat}.query.count\$" ],
        #   $nonhost_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
        addDEseq2( $config, $def, $summary_ref, "refseq_bacteria_${cat}_estimated", [ $refseq_bacteria_table, ".${cat}.estimated.count\$" ],
          $nonhost_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
        addDEseq2( $config, $def, $summary_ref, "refseq_bacteria_${cat}_aggregated", [ $refseq_bacteria_table, ".${cat}.aggregated.count\$" ],
          $nonhost_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
      }
      addDEseq2( $config, $def, $summary_ref, "refseq_bacteria_read", [ $refseq_bacteria_table, ".read.count\$" ],
        $nonhost_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
    }
  }

  my $nonhostXml   = [];
  my @overlap      = ();
  my @overlapNames = ();

  #Mapping unmapped reads to nonhost genome
  if ($search_nonhost_genome) {
    my $nonhost_genome_count_xml = [];
    my $microbial_genome_count_xml = [];
    for my $nonhostGroup (@nonhost_genome_groups) {
      addNonhostDatabase(
        $config, $def, $individual_ref, $summary_ref, "${nonhostGroup}_pm", $nonhost_genome_dir,    #general option
        $def->{"bowtie1_${nonhostGroup}_index"}, $identical_ref,                                    #bowtie option
        $def->{smallrnacount_option} . ' --keepChrInName --keepSequence',                                #count option
        $def->{nonhost_table_option} . ' --categoryMapFile ' . $def->{"${nonhostGroup}_species_map"},    #table option
        $identical_count_ref,
        $nonhostXml
      );

      addNonhostVis(
        $config, $def,
        $summary_ref,
        "nonhost_genome_${nonhostGroup}_vis",
        $data_visualization_dir,
        {
          rtemplate          => "countTableVisFunctions.R,countTableVis.R",
          output_file        => ".${nonhostGroup}Mapping.Result",
          output_file_ext    => ".Group.Piechart.png",
          parameterFile1_ref => [ "bowtie1_${nonhostGroup}_pm_table", ".Species.count\$" ],
          rCode              => 'maxCategory=4;' . $R_font_size,
        }
      );

      push @table_for_shortReadSource, ( "bowtie1_${nonhostGroup}_pm_table", ".read.count\$" );
      push @name_for_shortReadSource, ( $nonhostGroup);

      if ($do_comparison) {
        addDEseq2( $config, $def, $summary_ref, "${nonhostGroup}", [ "bowtie1_${nonhostGroup}_pm_table", ".Species.count\$" ],
          $nonhost_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
        addDEseq2( $config, $def, $summary_ref, "${nonhostGroup}_reads", [ "bowtie1_${nonhostGroup}_pm_table", ".read.count\$" ],
          $nonhost_genome_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
      }
      push @name_for_mapPercentage,      ( "bowtie1_${nonhostGroup}_pm_count", ".count.mapped.xml\$" );
      push @files_for_annotate_unmapped, ( "bowtie1_${nonhostGroup}_pm_count", ".count.mapped.xml\$" );
      push( @names_for_annotate_unmapped, $nonhostGroup );

      push @table_for_correlation, ( "bowtie1_${nonhostGroup}_pm_table", ".Species.count\$" );
      if ( $def->{read_correlation} ) {
        push @table_for_correlation, ( "bowtie1_${nonhostGroup}_pm_table", ".read.count\$", );
      }
      push @table_for_countSum,    ( "bowtie1_${nonhostGroup}_pm_table", ".Species.count\$" );
      push @table_for_readSummary, ( "bowtie1_${nonhostGroup}_pm_table", ".read.count\$" );
      push @overlap,               ( "bowtie1_${nonhostGroup}_pm_table", ".read.count\$" );
      
      push @$nonhost_count,( "bowtie1_${nonhostGroup}_pm_table", ".read.count\$" ); 

      push @$nonhost_genome_read_count,( "bowtie1_${nonhostGroup}_pm_table", ".read.count\$" ); 
      push @$nonhost_genome_count_xml, ( "bowtie1_${nonhostGroup}_pm_table", ".count.xml\$" );
      if( $nonhostGroup =~ /bacteria/ ){
        push @$bacteria_read_count, ( "bowtie1_${nonhostGroup}_pm_table", ".read.count\$" );
        push @$microbial_genome_count_xml, ( "bowtie1_${nonhostGroup}_pm_table", ".count.xml\$" );
      }

      my $nonhost_count2bam = ${nonhostGroup} . "_count2bam";
      if ( defined $def->{$nonhost_count2bam} && $def->{$nonhost_count2bam} ) {
        $config->{ "bowtie1_" . $nonhost_count2bam } = {
          class       => "CQS::Perl",
          perform     => 1,
          target_dir  => $nonhost_genome_dir . "/bowtie1_${nonhost_count2bam}",
          option      => getValue( $def, "${nonhost_count2bam}_option" ),
          output_ext  => getValue( $def, "${nonhost_count2bam}_output_ext" ),
          perlFile    => "bamByCount.pl",
          source_ref  => [ "bowtie1_${nonhostGroup}_pm", ".bam" ],
          source2_ref => [ "identical", ".dupcount\$" ],
          sh_direct   => 1,
          pbs         => {
            "nodes"    => "1:ppn=1",
            "walltime" => "2",
            "mem"      => "20gb"
          },
        };
        push( @$individual_ref, "bowtie1_" . $nonhost_count2bam );
      }
    }

    #print("perform_nonhost_genome_count=" . $perform_nonhost_genome_count . "\n");

    if ($perform_nonhost_genome_count){
      $config->{nonhost_genome_count} = {
        class              => "CQS::UniqueR",
        perform            => 1,
        target_dir         => $nonhost_genome_dir . "/nonhost_genome_count",
        option             => "",
        rtemplate          => "../SmallRNA/nonhostSampleCount.r",
        parameterSampleFile1_ref => $nonhost_genome_read_count,
        parameterSampleFile2 => {
          extension => ".nonhost_genome.tsv"
        },
        output_file_ext    => ".nonhost_genome.tsv",
        can_result_be_empty_file => 1,
        sh_direct          => 1,
        pbs                => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary_ref, "nonhost_genome_count" );

      $config->{microbial_genome_count} = {
        class              => "CQS::UniqueR",
        perform            => 1,
        target_dir         => $nonhost_genome_dir . "/microbial_genome_count",
        option             => "",
        rtemplate          => "../SmallRNA/nonhostSampleCount.r",
        parameterSampleFile1_ref => $bacteria_read_count,
        parameterSampleFile2 => {
          extension => ".microbial.tsv"
        },
        output_file_ext    => ".microbial.tsv",
        can_result_be_empty_file => 1,
        sh_direct          => 1,
        pbs                => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary_ref, "microbial_genome_count" );
    }
    
    push @name_for_readSummary, @nonhost_genome_group_names;
    push @overlapNames, @nonhost_genome_group_names;

    if ( $do_comparison && ( !$search_nonhost_genome_custom_group_only ) ) {
      addDeseq2Visualization( $config, $def, $summary_ref, "nonhost_genome",       \@nonhost_genome_groups,      $data_visualization_dir, "pairs_nonHostGroups_deseq2_vis_layout", $libraryKey );
      addDeseq2Visualization( $config, $def, $summary_ref, "nonhost_genome_reads", \@nonhost_genome_group_reads, $data_visualization_dir, "pairs_nonHostGroups_deseq2_vis_layout", $libraryKey );
    }
  }

  #Mapping unmapped reads to nonhost library
  if ($search_nonhost_library) {
    #Mapping unmapped reads to miRBase library
    addNonhostDatabase(
      $config, $def, $individual_ref, $summary_ref, "miRBase_pm", $nonhost_library_dir,    #general option
      $def->{bowtie1_miRBase_index}, $identical_ref,                                       #bowtie option
      $def->{mirbase_count_option} . " -m --keepChrInName --keepSequence",                 #count option
      $def->{nonhost_table_option},                                                        #table option
      $identical_count_ref,
      $nonhostXml
    );
    $config->{bowtie1_miRBase_pm_count}{can_result_be_empty_file} = 1;

    push @table_for_countSum, ( "bowtie1_miRBase_pm_table", "^((?!read|contig).)*\.count\$" );
    push @$nonhost_count, ( "bowtie1_miRBase_pm_table", '.read.count$' );

    #Mapping unmapped reads to tRNA library
    addNonhostDatabase(
      $config, $def, $individual_ref, $summary_ref, "tRNA_pm", $nonhost_library_dir,       #general option
      $def->{bowtie1_tRNA_index}, $identical_ref,                                          #bowtie option
      $def->{smallrnacount_option} . " --keepChrInName --keepSequence",                    #count option
      $def->{nonhost_table_option} . ' --categoryMapFile ' . $def->{trna_category_map},    #table option
      $identical_count_ref,
      $nonhostXml
    );
    push @$nonhost_count, ( "bowtie1_tRNA_pm_table", '.read.count$' );

    addNonhostVis(
      $config, $def,
      $summary_ref,
      "nonhost_library_tRNA_vis",
      $data_visualization_dir,
      {
        rtemplate          => "countTableVisFunctions.R,../SmallRNA/bacteriaTrnaMappingVis2.R",
        output_file        => ".tRNAMapping",
        output_file_ext    => ".Species12.csv;.tRNAType1.csv;.tRNAType2.csv",
        parameterFile1_ref => [ "bowtie1_tRNA_pm_table", ".count\$" ],
        parameterFile2     => getValue($def, "trna_category_map"),
        rCode              => 'maxCategory=3;' . $R_font_size,
      }
    );

    my $tRNA_bacteria_reads_task = "nonhost_tRNA_bacteria_reads";
    if ( $perform_nonhost_tRNA_bacteria_reads ) {
      my $folder            = $nonhost_library_dir . "/" . $tRNA_bacteria_reads_task;
      $config->{$tRNA_bacteria_reads_task} = {
        class              => "CQS::ProgramWrapper",
        perform            => 1,
        target_dir         => $folder,
        option             => "-c Bacteria",
        interpretor        => "python3",
        program            => "../SmallRNA/nonhostXmlCategoryRead.py",
        parameterFile1_arg => "-i",
        parameterFile1_ref => [ "bowtie1_tRNA_pm_table", ".Species.count" ],
        output_arg         => "-o",
        output_ext         => ".tRNA.Bacteria.read.count",
        sh_direct          => 1,
        pbs                => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary_ref, $tRNA_bacteria_reads_task );
    }

    if ( getValue( $def, "perform_nonhost_tRNA_coverage", 0 ) ) {
      my $positionTask      = "nonhost_library_tRNA_position";
      my $visualizationTask = $positionTask . "_vis_anticodon";
      my $folder            = $data_visualization_dir . "/" . $visualizationTask;
      $config->{$positionTask} = {
        class              => "CQS::ProgramWrapper",
        perform            => 1,
        target_dir         => $folder,
        option             => "-f " . $def->{bowtie1_tRNA_index} . ".fa -m " . getValue($def, "trna_map") . " -s " . $def->{nonhost_tRNA_coverage_species},
        interpretor        => "python3",
        program            => "../SmallRNA/tRNALibraryCoverage.py",
        parameterFile1_arg => "-i",
        parameterFile1_ref => [ "bowtie1_tRNA_pm_table", ".xml" ],
        output_arg         => "-o",
        output_ext         => ".tRNAlib.position",
        sh_direct          => 1,
        pbs                => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary_ref, $positionTask );

      addPositionVis(
        $config, $def,
        $summary_ref,
        $visualizationTask,
        $data_visualization_dir,
        {
          output_file        => ".nonhost_tRNAAnticodonPositionVis",
          parameterFile1_ref => [ $positionTask, ".position\$" ],
          rCode              => "countName<-\"TotalReads\""
        }
      );
    }

    #Mapping unmapped reads to rRNA library
    addNonhostDatabase(
      $config, $def, $individual_ref, $summary_ref, "rRNA_pm", $nonhost_library_dir,    #general option
      $def->{bowtie1_rRNA_index}, $identical_ref,                                       #bowtie option
      $def->{smallrnacount_option} . ' --keepChrInName --keepSequence --categoryMapFile ' . $def->{rrna_category_map},    #count option                                          #count option
      $def->{nonhost_table_option},                                                                                       #table option
      $identical_count_ref,
      $nonhostXml
    );
    push @$nonhost_count, ( "bowtie1_rRNA_pm_table", '.read.count$' );

    if ( getValue( $def, "perform_nonhost_rRNA_coverage" ) ) {
      my $positionTask      = "nonhost_library_rRNA_position";
      my $visualizationTask = $positionTask . "_vis";
      my $folder            = $data_visualization_dir . "/" . $visualizationTask;
      $config->{$positionTask} = {
        class                    => "CQS::ProgramWrapper",
        perform                  => 1,
        target_dir               => $folder,
        option                   => "-s " . getValue( $def, "nonhost_rRNA_coverage_species" ),
        interpretor              => "python3",
        program                  => "../SmallRNA/rRNALibraryCoverage.py",
        parameterSampleFile1_arg => "-i",
        parameterSampleFile1_ref => [ "bowtie1_rRNA_pm", ".bam" ],
        parameterSampleFile2_arg => "-c",
        parameterSampleFile2_ref => $identical_count_ref,
        output_arg               => "-o",
        output_ext               => ".rRNAlib.position",
        sh_direct                => 1,
        pbs                      => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push( @$summary_ref, $positionTask );

      addPositionVis(
        $config, $def,
        $summary_ref,
        $visualizationTask,
        $data_visualization_dir,
        {
          output_file        => ".nonhost_rRNAPositionVis",
          parameterFile1_ref => [ $positionTask, ".position\$" ],
          rCode              => "countName<-\"TotalReads\""
        }
      );
    }

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
        rCode              => 'maxCategory=NA;' . $R_font_size,
      }
    );

    push @table_for_shortReadSource, ( "bowtie1_tRNA_pm_table", ".read.count\$",  "bowtie1_rRNA_pm_table", ".read.count\$" );
    push @name_for_shortReadSource, ( "nonhost tRNA", "nonhost rRNA" );

    push( @name_for_mapPercentage,      "bowtie1_tRNA_pm_count", ".count.mapped.xml\$", "bowtie1_rRNA_pm_count", ".count.mapped.xml\$", );
    push( @files_for_annotate_unmapped, "bowtie1_tRNA_pm_count", ".count.mapped.xml\$", "bowtie1_rRNA_pm_count", ".count.mapped.xml\$", );
    push( @names_for_annotate_unmapped, "tRNA",                  "rRNA" );

    push @table_for_correlation, ( "bowtie1_tRNA_pm_table", $notReadCountPattern, "bowtie1_rRNA_pm_table", $notReadCountPattern );
    if ( $def->{read_correlation} ) {
      push @table_for_correlation, ( "bowtie1_tRNA_pm_table", ".read.count\$" );
      push @table_for_correlation, ( "bowtie1_rRNA_pm_table", ".read.count\$" );
      if($perform_nonhost_tRNA_bacteria_reads){
        push @table_for_correlation, ( $tRNA_bacteria_reads_task, ".read.count\$" );
      }
    }
    push @table_for_countSum,    ( "bowtie1_tRNA_pm_table", ".Species.count\$", "bowtie1_rRNA_pm_table", "$task_name\.count\$" );
    push @table_for_readSummary, ( "bowtie1_tRNA_pm_table", ".read.count\$",     "bowtie1_rRNA_pm_table", ".read.count\$" );
    push @name_for_readSummary,  ( "Non host tRNA",         "Non host rRNA" );
    push @overlap,               ( "bowtie1_tRNA_pm_table", ".read.count\$",     "bowtie1_rRNA_pm_table", ".read.count\$" );
    push @overlapNames,          ( "Non host tRNA",         "Non host rRNA" );
    if ($do_comparison) {
      my $tRNADeseq2 = [];
      push @$tRNADeseq2,
        addDEseq2( $config, $def, $summary_ref, "nonhost_tRNA", [ "bowtie1_tRNA_pm_table", ".count\$" ], $nonhost_library_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
      push @$tRNADeseq2,
        addDEseq2( $config, $def, $summary_ref, "nonhost_tRNA_species", [ "nonhost_library_tRNA_vis", ".Species12.csv\$" ],
        $nonhost_library_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
      push @$tRNADeseq2,
        addDEseq2( $config, $def, $summary_ref, "nonhost_tRNA_type", [ "nonhost_library_tRNA_vis", ".tRNAType1.csv\$" ], $nonhost_library_dir, $DE_min_median_read_smallRNA, $libraryFile,
        $libraryKey );
      push @$tRNADeseq2,
        addDEseq2( $config, $def, $summary_ref, "nonhost_tRNA_anticodon", [ "nonhost_library_tRNA_vis", ".tRNAType2.csv\$" ],
        $nonhost_library_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
      push @$tRNADeseq2,
        addDEseq2( $config, $def, $summary_ref, "nonhost_tRNA_reads", [ "bowtie1_tRNA_pm_table", ".read.count\$" ], $nonhost_library_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
      addDEseq2( $config, $def, $summary_ref, "nonhost_tRNA_category", [ "bowtie1_tRNA_pm_table", ".Species.count\$" ], $nonhost_library_dir, $DE_min_median_read_smallRNA, $libraryFile,
        $libraryKey );

      if($perform_nonhost_tRNA_bacteria_reads){
        addDEseq2( $config, $def, $summary_ref, $tRNA_bacteria_reads_task, [ $tRNA_bacteria_reads_task, ".read.count\$" ], $nonhost_library_dir, $DE_min_median_read_smallRNA, $libraryFile,
          $libraryKey );
      }


      addDeseq2Visualization( $config, $def, $summary_ref, "nonhost_tRNA", [ "nonhost_tRNA", "nonhost_tRNA_species", "nonhost_tRNA_type", "nonhost_tRNA_anticodon", "nonhost_tRNA_reads" ],
        $data_visualization_dir, "pairs_nonHostLibrary_deseq2_vis_layout", $libraryKey );

      addDEseq2( $config, $def, $summary_ref, "nonhost_rRNA", [ "bowtie1_rRNA_pm_table", ".count\$" ], $nonhost_library_dir, $DE_min_median_read_smallRNA, $libraryFile, $libraryKey );
    }
  }

  my $perform_nonhost_count = ($search_nonhost_genome || $search_nonhost_library) && getValue($def, "perform_nonhost_count", 1);
  if ($perform_nonhost_count) {
    $config->{nonhost_count} = {
      class              => "CQS::UniqueR",
      perform            => 1,
      target_dir         => $data_visualization_dir . "/nonhost_count",
      option             => "",
      rtemplate          => "../SmallRNA/nonhostSampleCount.r",
      parameterSampleFile1_ref => $nonhost_count,
      parameterSampleFile2 => {
        extension => ".nonhost.tsv"
      },
      output_file_ext    => ".nonhost.tsv",
      can_result_be_empty_file => 1,
      sh_direct          => 1,
      pbs                => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    };
    push( @$summary_ref, "nonhost_count" );
  }

  if ( $def->{perform_nonhost_mappedToHost} ) {
    my $bowtie1readTask = "bowtie1_nonhost_mappedreads";
    $config->{$bowtie1readTask} = {
      class                    => "CQS::ProgramIndividualWrapper",
      perform                  => 1,
      target_dir               => "${data_visualization_dir}/$bowtie1readTask",
      option                   => "",
      interpretor              => "python3",
      program                  => "../SmallRNA/nonhostXmlToFastq.py",
      source_arg               => "-i",
      source_ref               => $nonhostXml,
      parameterSampleFile2_arg => "-f",
      parameterSampleFile2_ref => $identical_ref,
      output_to_same_folder    => 1,
      output_arg               => "-o",
      output_ext               => ".fastq.gz",
      join_arg                 => 1,
      sh_direct                => 1,
      pbs                      => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "10gb"
      },
    };
    push @$individual_ref, $bowtie1readTask;

    my $bowtie1readMapTask = "bowtie1_nonhost_mappedreads_host";
    addBowtie( $config, $def, $individual_ref, $bowtie1readMapTask, $data_visualization_dir, $def->{bowtie1_index}, [$bowtie1readTask], $def->{bowtie1_option_2mm} );

    my $bowtie1readMapMismatchTask = "bowtie1_nonhost_mappedreads_host_mismatch_table";
    $config->{$bowtie1readMapMismatchTask} = {
      class                    => "CQS::ProgramWrapper",
      perform                  => 1,
      target_dir               => "${data_visualization_dir}/$bowtie1readMapMismatchTask",
      option                   => "-m 2",
      interpretor              => "python3",
      program                  => "../SmallRNA/bamMismatchTable.py",
      parameterSampleFile1_arg => "-i",
      parameterSampleFile1_ref => [ $bowtie1readMapTask, ".bam\$" ],
      parameterSampleFile2_arg => "-f",
      parameterSampleFile2_ref => [ $bowtie1readTask, ".fastq.gz\$" ],
      parameterSampleFile3_arg => "-c",
      parameterSampleFile3_ref => $identical_count_ref,
      output_to_same_folder    => 1,
      output_arg               => "-o",
      output_ext               => ".tsv",
      sh_direct                => 1,
      pbs                      => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "10gb"
      },
    };
    push @$summary_ref, $bowtie1readMapMismatchTask;
  }

  if ( $search_nonhost_database && ( !$search_nonhost_genome_custom_group_only ) ) {
    if ($perform_nonhost_overlap_vis) {
      
      my $rcode = "categoriesNames=c(\"" . join("\",\"", @overlapNames) . "\");";
      
      $config->{nonhost_overlap_vis} = {
        class                     => "CQS::UniqueR",
        perform                   => 1,
        target_dir                => $data_visualization_dir . "/nonhost_overlap_vis",
        rtemplate                 => "countTableVisFunctions.R,NonHostOverlap.R",
        output_file               => ".NonHost.Reads",
        output_file_ext           => ".Overlap.csv;.Barplot.png;",
        parameterSampleFile1_ref  => \@overlap,
        parameterSampleFile2Order => $def->{groups_order},
        parameterSampleFile2      => $groups,
        parameterSampleFile3      => $def->{groups_vis_layout},
        parameterFile3_ref        => [ "fastqc_count_vis", ".Reads.csv\$" ],
        sh_direct                 => 1,
        rCode                     => 'maxCategory=8;' . $rcode . $R_font_size,
        pbs                       => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push @$summary_ref, "nonhost_overlap_vis";
    }

    $config->{final_unmapped_reads} = {
      class            => "CQS::ProgramWrapperOneToOne",
      perform          => 1,
      target_dir       => $nonhost_blast_dir . "/final_unmapped_reads",
      option           => "-o __NAME__.unmapped.fastq.gz",
      interpretor      => "python3",
      program          => "../SmallRNA/unmappedReadsToFastq.py",
      source_arg       => "-i",
      source_ref       => $identical_ref,
      parameterSampleFile2_arg      => "-m",
      parameterSampleFile2_ref      => $nonhost_count,
      parameterSampleFile3_arg      => "-c",
      parameterSampleFile3_ref      => $identical_count_ref,
      output_to_same_folder => 1,
      output_arg       => "-o",
      output_file_ext  => ".unmapped.fastq.gz,.unmapped.fastq.dupcount,.unmapped.fastq.gz.info",
      sh_direct        => 1,
      pbs              => {
        "nodes"    => "1:ppn=1",
        "walltime" => "1",
        "mem"      => "10gb"
      },
    };
    $config->{final_unmapped_reads_summary} = {
      class      => "CQS::CQSDatatable",
      perform    => 1,
      target_dir => $nonhost_blast_dir . "/final_unmapped_reads_summary",
      source_ref => [ "final_unmapped_reads", ".unmapped.fastq.gz.info\$" ],
      option     => "",
      sh_direct  => 1,
      pbs        => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "40gb"
      },
    };

    $identical_ref = [ "final_unmapped_reads", ".fastq.gz\$" ];
    push @$individual_ref,      ("final_unmapped_reads");
    push @$summary_ref,         ("final_unmapped_reads_summary");
    push @table_for_pieSummary, ( "final_unmapped_reads", ".dupcount" );
    push @name_for_pieSummary,  ("UnMapped");

    if ( $def->{perform_map_percentage} ) {
      $config->{map_percentage} = {
        class      => "SmallRNA::MapPercentage",
        perform    => 1,
        target_dir => $data_visualization_dir . "/map_percentage",
        source_ref => \@name_for_mapPercentage,
        sh_direct  => 0,
        pbs        => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "12",
          "mem"       => "10gb"
        },
      };
      push @$individual_ref, "map_percentage";
    }
  }

  if ($perform_annotate_unmapped_reads) {
    $config->{annotate_unmapped_reads} = {
      class            => "SmallRNA::AnnotateUnmappedReads",
      perform          => 1,
      target_dir       => $nonhost_blast_dir . "/annotate_unmapped_reads",
      source_ref       => \@reads_for_annoate_unmapped,
      mapped_files_ref => \@files_for_annotate_unmapped,
      mapped_names     => join( ',', @names_for_annotate_unmapped ),
      min_count        => 2,
      sh_direct        => 1,
      pbs              => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "12",
        "mem"       => "10gb"
      },
    };
    push @$individual_ref, "annotate_unmapped_reads";
  }

  my $corr_output_file_ext      = ".Correlation.png;.heatmap.png;.PCA.png;";
  my $corr_output_file_task_ext = "";
  if ( ( defined $def->{groups} ) && ( scalar( keys %$groups ) >= 3 ) ) {
    $corr_output_file_task_ext = ".Group.heatmap.png;.Group.Correlation.Cluster.png;";
  }

  my $correlation_folder = "count_table_correlation";
  if($normalize_by ne 'None'){
    $correlation_folder = $correlation_folder . '_' . $normalize_by;
  }
  $config->{count_table_correlation} = {
    class                     => "CQS::CountTableGroupCorrelation",
    perform                   => 1,
    target_dir                => "$data_visualization_dir/$correlation_folder",
    rtemplate                 => "countTableVisFunctions.R,countTableGroupCorrelation.R",
    output_file               => "parameterSampleFile1",
    output_file_ext           => $corr_output_file_ext,
    output_file_task_ext      => $corr_output_file_task_ext,
    output_to_result_dir      => getValue($def, "correlation_output_to_result_dir", 0),
    output_include_folder_name => getValue($def, "correlation_output_include_folder_name", 1),
    parameterSampleFile1_ref  => \@table_for_correlation,
    parameterSampleFile2      => $def->{tRNA_vis_group},
    parameterSampleFile2Order => $def->{groups_order},
    parameterFile3_ref        => [ "fastqc_count_vis", ".Reads.csv\$" ],
    rCode                     => $def->{correlation_rcode} . $R_font_size,
    sh_direct                 => 1,
    pbs                       => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "20gb"
    },
  };

  if($search_host_genome){
    $config->{count_table_correlation}{parameterFile2_ref} = [ "bowtie1_genome_1mm_NTA_smallRNA_category", ".Category.Table.csv" ];
  }

  if ( defined $def->{groups}) {
    if(defined $def->{correlation_groups}){
      $config->{count_table_correlation}{parameterSampleFile2} = $def->{correlation_groups};
    }elsif ( defined $def->{correlation_groups_dic} ) {
      my $correlationGroups = get_pair_group_sample_map( $def->{correlation_groups_dic}, $def->{groups} );
      if ( getValue( $def, "correlation_all", 1 ) and (not defined $correlationGroups->{all}) ) {
        $correlationGroups->{all} = $def->{groups};
      }
      $config->{count_table_correlation}{parameterSampleFile2} = $correlationGroups;
    }
  }

  push @$summary_ref, ("count_table_correlation");

  if(defined $deseq2Task){
    if($def->{perform_permanova}){
      my $log_transform = getValue($def, "permanova_log_transform", 1);
      my $log_prefix = $log_transform ? ".log2" : "";

      $config->{count_table_permanova} = {
        class                     => "CQS::UniqueR",
        perform                   => 1,
        target_dir                => $data_visualization_dir . "/count_table_permanova",
        rtemplate                 => "../SmallRNA/permanova.r",
        output_file               => "",
        output_file_ext           => "$log_prefix.PCoA.pdf",
        output_file_task_ext      => "$log_prefix.permanova.txt;$log_prefix.betadisper.txt",
        parameterSampleFile1_ref  => \@table_for_correlation,
        parameterSampleFile2      => $config->{count_table_correlation}{parameterSampleFile2},
        parameterSampleFile2Order => $def->{groups_order},
        parameterSampleFile3_ref  => [ $deseq2Task, ".design\$"],
        parameterSampleFile4  => {
          log_transform => getValue($def, "permanova_log_transform", 1)
        },
        parameterFile2_ref        => $config->{count_table_correlation}{parameterFile2_ref},
        parameterFile3_ref        => $config->{count_table_correlation}{parameterFile3_ref},
        rCode                     => $config->{count_table_correlation}{rCode},
        sh_direct                 => 1,
        pbs                       => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };

      push @$summary_ref, ("count_table_permanova");
    }
  }

  my $paramFile = undef;
  if(defined $config->{nonhost_genome_count}){
    $paramFile = ["nonhost_genome_count"];
  }

  if(scalar(@table_for_countSum) > 0) {
    $config->{reads_in_tasks} = {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      target_dir               => $data_visualization_dir . "/reads_in_tasks",
      rtemplate                => "countTableVisFunctions.R,ReadsInTasks.R",
      output_file_ext          => ".TaskReads.csv",
      parameterSampleFile1_ref => \@table_for_countSum,
      parameterFile3_ref       => [ "fastqc_count_vis", ".Reads.csv\$" ],
      rCode                    => $R_font_size,
      sh_direct                => 1,
      pbs                      => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "12",
        "mem"       => "10gb"
      },
    };

    push @$summary_ref, ("reads_in_tasks");
  }

  if ( $search_host_genome && $search_nonhost_database ) {
    my $cur_r_code = $R_font_size . " " . getUniqueGroups($def);
    $config->{reads_in_tasks_pie} = {
      class                => "CQS::UniqueR",
      suffix               => "_pie",
      perform              => 1,
      target_dir           => $data_visualization_dir . "/reads_in_tasks",
      rtemplate            => "countTableVisFunctions.R,ReadsInTasksPie.R",
      output_file_ext      => ".NonParallel.TaskReads.csv",
      parameterFile1_ref   => [ "bowtie1_genome_1mm_NTA_smallRNA_info", ".mapped.count\$" ],
      parameterFile2_ref   => [ "final_unmapped_reads_summary", ".count\$" ],
      parameterFile3_ref   => $paramFile,
      parameterSampleFile1 => $groups,
      parameterSampleFile2 => $def->{groups_vis_layout},
      rCode                => $cur_r_code,
      sh_direct            => 1,
      pbs                  => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "12",
        "mem"       => "10gb"
      },
    };
    $config->{reads_in_tasks_all} = {
      class              => "CQS::UniqueR",
      suffix             => "_all",
      perform            => 1,
      target_dir         => $data_visualization_dir . "/reads_in_tasks",
      rtemplate          => "countTableVisFunctions.R,ReadsInTasksAll.R",
      output_file_ext    => ".All.TaskReads.csv",
      parameterFile1_ref => [ "reads_in_tasks", ".TaskReads.csv\$" ],
      parameterFile2_ref => [ "reads_in_tasks_pie", ".NonParallel.TaskReads.csv\$" ],
      parameterFile3_ref => [ "fastqc_count_vis", ".Reads.csv\$" ],
      rCode              => $R_font_size,
      sh_direct          => 1,
      pbs                => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "12",
        "mem"       => "10gb"
      },
    };
    push @$summary_ref, ( "reads_in_tasks_pie", "reads_in_tasks_all" );

    if($perform_nonhost_genome_count){  
      my $rCode = getUniqueGroups($def);
      $config->{host_microbial_vis} = {
        class                     => "CQS::UniqueR",
        perform                   => 1,
        target_dir                => $data_visualization_dir . "/host_microbial_vis",
        rtemplate                 => "../SmallRNA/hostMicrobialVis.r",
        output_file               => ".reads",
        output_file_ext           => ".pdf",
        parameterSampleFile1      => $groups,
        parameterSampleFile2      => $def->{groups_vis_layout},
        parameterFile1_ref        => [ "reads_in_tasks_pie", ".NonParallel.TaskReads.csv"],
        parameterFile2_ref        => [ "microbial_genome_count", ".microbial.tsv\$" ],
        sh_direct                 => 1,
        rCode                     => $rCode,
        pbs                       => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      };
      push @$summary_ref, "host_microbial_vis";
    }
  }

  if ($perform_class_independent_analysis) {
    if(scalar(@table_for_readSummary) > 0){
      my $name_for_readSummary_r = "readFilesModule=c('" . join( "','", @name_for_readSummary ) . "'); ";
      $config->{top_sequence_mapped_in_categories} = {
        class                    => "CQS::UniqueR",
        perform                  => 1,
        target_dir               => $data_visualization_dir . "/top_sequence_mapped_in_categories",
        rtemplate                => "countTableVisFunctions.R,TopReadsMappingSummary.R",
        output_file_ext          => ".ReadsMapping.Summary.csv",
        parameterFile1_ref       => [ "identical_sequence_count_table", $task_name . "_sequence.read.count\$" ],
        parameterSampleFile1_ref => \@table_for_readSummary,
        parameterSampleFile2     => $groups,
        parameterSampleFile3     => $def->{groups_vis_layout},
        rCode                    => $name_for_readSummary_r . $R_font_size,
        sh_direct                => 1,
        pbs                      => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "12",
          "mem"       => "10gb"
        },
      };
      push @$summary_ref, "top_sequence_mapped_in_categories";
    }
    
    if ($perform_short_reads_source) {
      $config->{short_reads_source} = {
        'class'                    => 'CQS::ProgramWrapper',
        'parameterSampleFile1_arg' => '-m',
        'parameterSampleFile1_ref' => [ "bowtie1_genome_1mm_NTA", ".bam.max.txt" ],
        'parameterSampleFile2_arg' => '-a',
        'parameterSampleFile2_ref' => \@table_for_readSummary,
        'parameterFile1_arg'       => "-i",
        'parameterFile1_ref'       => [ $short_reads_table, ".count.txt\$" ],
        'option'                   => "-n \"" . join( ",", @name_for_readSummary ) . "\"",
        'interpretor'              => 'python3',
        'program'                  => '../SmallRNA/shortReadSource.py',
        'target_dir'               => $data_visualization_dir . "/short_reads_source",
        'output_ext'               => '.tsv',
        'output_arg'               => '-o',
        'output_to_same_folder'    => 1,
        'sh_direct'                => 1,
        'perform'                  => 1,
        'pbs'                      => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "2",
          "mem"       => "10gb"
        },
      };
      push @$summary_ref, "short_reads_source";
    }
  }

  #add time cost task in the end of pipeline
  #search not identical reads to genome, for IGV
  if ( $search_host_genome && $search_not_identical ) {
    addBowtie( $config, $def, $individual_ref, "bowtie1_genome_1mm_notidentical", $host_genome_dir, $def->{bowtie1_index}, $not_identical_ref, $def->{bowtie1_option_1mm} );
  }

  #blast top reads
  if ( $blast_top_reads and $perform_class_independent_analysis ) {
    my $deseq2TopTask = getDEseq2TaskName( "top${top_read_number}_minicontigs", $libraryKey, $def );
    if ($do_comparison) {
      addDeseq2SignificantSequenceBlastn( $config, $def, $summary_ref, $deseq2TopTask, $class_independent_dir );
    }
    else {
      #addBlastn( $config, $def, $summary_ref, "identical_sequence_top${top_read_number}_contig_blast",     "identical_sequence_count_table", "sequence.count.fasta\$",   $class_independent_dir );
      #addBlastn( $config, $def, $summary_ref, "identical_sequence_top${top_read_number}_read_blast",       "identical_sequence_count_table", "read.count.fasta\$",       $class_independent_dir );
      addBlastn( $config, $def, $summary_ref, "identical_sequence_top${top_read_number}_minicontig_blast", "identical_sequence_count_table", "minicontig.count.fasta\$", $class_independent_dir );
    }
  }

  #blast unmapped reads
  if ($blast_unmapped_reads) {
    $config->{"final_unmapped_reads_table"} = {
      class         => "CQS::SmallRNASequenceCountTable",
      perform       => 1,
      target_dir    => $nonhost_blast_dir . "/final_unmapped_reads_table",
      option        => getValue( $def, "sequence_count_option" ),
      source_ref    => [ $identical_ref->[0], ".dupcount\$" ],
      suffix        => "_unmapped",
      sh_direct     => 1,
      cluster       => $cluster,
      groups        => $def->{groups},
      deseq2_groups => $def->{deseq2_groups},
      pairs         => $def->{pairs},
      pbs           => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "10gb"
      },
    };
    push @$summary_ref, "final_unmapped_reads_table";

    if ($do_comparison) {
      $deseq2Task = addDEseq2(
        $config, $def, $summary_ref,
        "final_unmapped_reads_minicontigs",
        [ "final_unmapped_reads_table", ".minicontig.count\$" ],
        $nonhost_blast_dir, $DE_min_median_read_top, $libraryFile, $libraryKey
      );
      addDeseq2SignificantSequenceBlastn( $config, $def, $summary_ref, $deseq2Task, $nonhost_blast_dir );
    }
    else {
      addBlastn( $config, $def, $summary_ref, "final_unmapped_reads_minicontigs_blast", "final_unmapped_reads_table", "minicontig.count.fasta\$", $nonhost_blast_dir );
    }
  }

  #tDRmapper
  if ($perform_tDRmapper) {
    my $tools_dir = create_directory_or_die( $def->{target_dir} . "/other_tools" );
    $config->{"tDRmapper"} = {
      class      => "CQS::Perl",
      perform    => 1,
      target_dir => $tools_dir . "/tDRmapper",
      perlFile   => "runtDRmapper.pl",
      option     => $def->{tDRmapper} . " " . $def->{tDRmapper_fasta},
      source_ref => $not_identical_ref,
      output_ext => "_clipped_identical.fastq.hq_cs",
      sh_direct  => 0,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "40gb"
      },
    };
    push @$individual_ref, ("tDRmapper");
  }

  #check batch effect
  if ( defined $batchGroups ) {
    for my $batchGroup ( sort keys %$batchGroups ) {
      my $batchName = "count_table_correlation" . "_" . $batchGroup;

      my $batchConfig = dclone( $config->{"count_table_correlation"} );
      $batchConfig->{target_dir}                = $batch_dir . "/" . $batchName;
      $batchConfig->{suffix}                    = "_" . $batchGroup;
      $batchConfig->{output_to_result_dir}      = "1";
      $batchConfig->{parameterSampleFile2}      = $batchGroups->{$batchGroup};
      $batchConfig->{parameterSampleFile2Order} = undef;
      $batchConfig->{rCode}                     = ( defined $batchConfig->{rCode} ? $batchConfig->{rCode} : "" ) . "visLayoutAlphabet=TRUE;" . $R_font_size;
      $config->{$batchName}                     = $batchConfig;

      push @$summary_ref, ($batchName);
    }

    #fastq_len_vis
    for my $batchGroup ( sort keys %$batchGroups ) {
      my $batchName = "fastq_len_vis" . "_" . $batchGroup;

      my $batchConfig = dclone( $config->{"fastq_len_vis"} );
      $batchConfig->{target_dir}           = $batch_dir . "/" . $batchName;
      $batchConfig->{parameterSampleFile2} = $batchGroups->{$batchGroup};
      $batchConfig->{parameterSampleFile3} = $batchLayout->{$batchGroup};
      $batchConfig->{output_file}          = ".len_" . $batchGroup;
      $batchConfig->{rCode}                = ( defined $batchConfig->{rCode} ? $batchConfig->{rCode} : "" ) . "visLayoutAlphabet=TRUE;" . $R_font_size;

      $config->{$batchName} = $batchConfig;

      push @$summary_ref, ($batchName);
    }

    #reads_in_tasks_pie
    for my $batchGroup ( sort keys %$batchGroups ) {
      my $batchName = "reads_in_tasks_pie" . "_" . $batchGroup;

      my $batchConfig = dclone( $config->{"reads_in_tasks_pie"} );
      $batchConfig->{target_dir}                 = $batch_dir . "/" . $batchName;
      $batchConfig->{parameterSampleFile2}       = $batchGroups->{$batchGroup};
      $batchConfig->{parameterSampleFile2_order} = undef;
      $batchConfig->{parameterSampleFile3}       = $batchLayout->{$batchGroup};
      $batchConfig->{output_file}                = ".reads_" . $batchGroup;
      $batchConfig->{rCode}                      = ( defined $batchConfig->{rCode} ? $batchConfig->{rCode} : "" ) . "visLayoutAlphabet=TRUE;";

      $config->{$batchName} = $batchConfig;

      push @$summary_ref, ($batchName);
    }

    #bowtie1_genome_1mm_NTA_smallRNA_category
    for my $batchGroup ( sort keys %$batchGroups ) {
      my $batchName = "bowtie1_genome_1mm_NTA_smallRNA_category" . "_" . $batchGroup;

      my $batchConfig = dclone( $config->{"bowtie1_genome_1mm_NTA_smallRNA_category"} );
      $batchConfig->{target_dir}                = $batch_dir . "/" . $batchName;
      $batchConfig->{parameterSampleFile2}      = $batchGroups->{$batchGroup};
      $batchConfig->{parameterSampleFile2Order} = undef;
      $batchConfig->{parameterSampleFile3}      = $batchLayout->{$batchGroup};
      $batchConfig->{rCode}                     = ( defined $batchConfig->{rCode} ? $batchConfig->{rCode} : "" ) . "drawInvidividual=FALSE;visLayoutAlphabet=TRUE;";

      $config->{$batchName} = $batchConfig;

      push @$summary_ref, ($batchName);
    }
  }

  if ( $config->{fastqc_count_vis} && $config->{bowtie1_genome_1mm_NTA_smallRNA_category} ) {

    my $read_in_tasks_file = undef;
    if( $config->{reads_in_tasks_pie}){
      $read_in_tasks_file = [ "reads_in_tasks_pie", ".NonParallel.TaskReads.csv\$" ];
    }
    $config->{read_summary} = {
      class              => "CQS::UniqueR",
      perform            => 1,
      target_dir         => $data_visualization_dir . "/read_summary",
      rtemplate          => "../SmallRNA/readSummary.R",
      #output_file_ext    => ".perc.png;.count.png",
      output_file_ext    => ".png",
      parameterFile1_ref => [ "fastqc_count_vis", ".countInFastQcVis.Result.Reads.csv\$" ],
      parameterFile2_ref => $read_in_tasks_file,
      parameterFile3_ref => [ "bowtie1_genome_1mm_NTA_smallRNA_category", ".Category.Table.csv\$" ],
      rCode              => "",
      sh_direct          => 1,
      pbs                => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "12",
        "mem"       => "10gb"
      },
    };
    push @$summary_ref, ("read_summary");
  }

  if ( $search_host_genome && $search_nonhost_database && $perform_bacteria_count ) {
    $config->{bacteria_count} = {
      'class'                    => 'CQS::ProgramWrapper',
      'parameterSampleFile1_arg' => '-g',
      'parameterSampleFile1_ref' => $bacteria_read_count,
      'parameterFile2_arg' => '-d',
      'parameterFile2_ref' => [ "bowtie1_rRNA_pm_table", ".count.xml\$" ],
      'parameterFile3_arg'       => '-t',
      'parameterFile3_ref'       => [ "reads_in_tasks_pie", ".NonParallel.TaskReads.csv\$" ],
      'option'                   => "",
      'interpretor'              => 'python3',
      'program'                  => '../SmallRNA/getBacteriaCount.py',
      'target_dir'               => $data_visualization_dir . "/bacteria_count",
      'output_file_ext'          => '.tsv;.tsv.summary;.tsv.summary.png;.tsv.summary.rpm.csv',
      'output_arg'               => '-o',
      'output_to_same_folder'    => 1,
      'sh_direct'                => 1,
      'perform'                  => 1,
      'pbs'                      => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "2",
        "mem"       => "10gb"
      },
    };
    push @$summary_ref, "bacteria_count";

    if(getValue($def, "search_refseq_bacteria")){
      $config->{bacteria_count_summary} = {
        'class'                    => 'CQS::UniqueR',
        'parameterFile1_ref' => [ "bowtie1_bacteria_group1_pm_table", ".read.count\$" ],
        'parameterFile2_ref' => [ "bowtie1_bacteria_group2_pm_table", ".read.count\$" ],
        'parameterFile3_ref' => [ "refseq_bacteria_table", ".read.count\$" ],
        'option'                   => "",
        'rtemplate'                => 'countTableVisFunctions.R,../SmallRNA/bacteriaReadSummary.R',
        'target_dir'               => $data_visualization_dir . "/bacteria_count_summary",
        'output_file_ext'          => '.csv,.bar1.png,.bar2.png',
        'sh_direct'                => 1,
        'perform'                  => 1,
        'pbs'                      => {
          "email"     => $def->{email},
          "emailType" => $def->{emailType},
          "nodes"     => "1:ppn=1",
          "walltime"  => "2",
          "mem"       => "10gb"
        },
      };
      push @$summary_ref, "bacteria_count_summary";
    }
  }

  if (getValue($def, "perform_search_fasta", 0)){
    my $old_use_intermediate_dir = $def->{"use_intermediate_dir"};

    $def->{"use_intermediate_dir"} = 0;

    my $fasta = getValue($def, "search_fasta_file");
    my $fasta_name = getValue($def, "search_fasta_name");
    my $search_fasta_folder = getValue($def, "search_fasta_folder", $def->{target_dir} . "/$fasta_name");
    my $fasta_folder = create_directory_or_die($search_fasta_folder);

    #my $bowtie1_option_1mm = "-a -m 100 --best --strata -v 1 -l 12 -p 8";
    my $bowtie1_option_pm = "-a -m 100 --best --strata -v 0 -l 12 -p 8";

    my $bowtie_index_task;

    my $search_fasta_tasks = [];

    if ( -e $fasta . ".4.ebwt") {
      $bowtie_index_task = $fasta;
    }else{
      $bowtie_index_task = "bowtie1_" . $fasta_name . "_00_index";
      add_bowtie_index($config, $def, $individual_ref, $fasta_folder, $bowtie_index_task, $fasta);
      push(@$search_fasta_tasks, $bowtie_index_task);
    }

    if (getValue($def, "search_fasta_all_reads", 0)){
      my @tasks = add_search_fasta($config, $def, $individual_ref, $summary_ref, $fasta_name . "_all", $fasta_folder, $bowtie_index_task, ["identical", "fastq.gz"], $bowtie1_option_pm, "", "" );
      push(@$search_fasta_tasks, @tasks)
    }

    if (getValue($def, "search_fasta_nonhost_reads", 0)){
      my @tasks =add_search_fasta($config, $def, $individual_ref, $summary_ref, $fasta_name . "_nonhost", $fasta_folder, $bowtie_index_task, ["bowtie1_genome_unmapped_reads", '.unmapped.fastq.gz$'], $bowtie1_option_pm, "", "" );
      push(@$search_fasta_tasks, @tasks)
    }

    if(getValue($def, "search_fasta_sequencetask", 1)){
      $config->{search_fasta_sequencetask} = {
        class      => "CQS::SequenceTaskNoDependence",
        perform    => 1,
        target_dir => $fasta_folder . "/search_fasta_sequencetask",
        option     => "",
        source     => {
          step1 => $search_fasta_tasks,
        },
        sh_direct => 1,
        pbs       => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "10",
          "mem"       => "10gb"
        },
      };
    }

    $def->{"use_intermediate_dir"} = $old_use_intermediate_dir;
  }

  if ( getValue( $def, "perform_report" ) ) {
    my @report_files = ();
    my @report_names = ();
    my @copy_files   = ();

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

    if ( defined $config->{read_summary} ) {
      push( @report_files, "read_summary", ".png" );
      push( @report_names, "read_summary" );
    }

    if ( defined $config->{fastq_len} ) {
      push( @report_files, "fastq_len_vis", ".lengthDistribution.png" );
      push( @report_names, "fastq_len" );
    }

    if ( defined $config->{bowtie1_genome_1mm_NTA_smallRNA_category} ) {
      if ( !defined $config->{read_summary} ) {
        push( @report_files, "bowtie1_genome_1mm_NTA_smallRNA_category", ".Category1.Barplot.png" );
        push( @report_files, "bowtie1_genome_1mm_NTA_smallRNA_category", ".Category2.Barplot.png" );
        push( @report_names, "category_mapped_bar",                      "category_smallrna_bar" );
      }

      if (defined $def->{groups}){
        push( @report_files, "bowtie1_genome_1mm_NTA_smallRNA_category", ".Category1.Group.Piechart.png" );
        push( @report_files, "bowtie1_genome_1mm_NTA_smallRNA_category", ".Category2.Group.Piechart.png" );
        push( @report_names, "category_mapped_group",                    "category_smallrna_group" );
      }
    }

    if ( defined $config->{count_table_correlation} ) {
      my $hasGroupHeatmap = 0;
      if ( defined $def->{groups} ) {
        my $groups = $def->{groups};
        $hasGroupHeatmap = scalar( keys %$groups ) > 2;
      }
      if ( defined $config->{bowtie1_genome_1mm_NTA_smallRNA_table} ) {
        push( @report_files, "count_table_correlation",   ".miRNA.count.heatmap.png" );
        push( @report_files, "count_table_correlation",   ".miRNA.count.PCA.png" );
        push( @report_names, "correlation_mirna_heatmap", "correlation_mirna_pca" );

        if ($hasGroupHeatmap) {
          push( @report_files, "count_table_correlation",         ".miRNA.count.Group.heatmap.png" );
          push( @report_files, "count_table_correlation",         ".miRNA.count.Group.Correlation.Cluster.png" );
          push( @report_names, "correlation_mirna_group_heatmap", "correlation_mirna_corr_cluster" );
        }

        if(not $hasMicroRNAOnly){
          push( @report_files, "count_table_correlation",  'smallRNA_1mm_.+.tRNA.count.heatmap.png' );
          push( @report_files, "count_table_correlation",  'smallRNA_1mm_.+.tRNA.count.PCA.png' );
          push( @report_names, "correlation_trna_heatmap", "correlation_trna_pca" );

          push( @report_files, "count_table_correlation",  "smallRNA_1mm_.+.rRNA.count.heatmap.png" );
          push( @report_files, "count_table_correlation",  "smallRNA_1mm_.+.rRNA.count.PCA.png" );
          push( @report_names, "correlation_rrna_heatmap", "correlation_rrna_pca" );

          if($def->{hasYRNA}){
            push( @report_files, "count_table_correlation",  "smallRNA_1mm_.+.yRNA.count.heatmap.png" );
            push( @report_files, "count_table_correlation",  "smallRNA_1mm_.+.yRNA.count.PCA.png" );
            push( @report_names, "correlation_yrna_heatmap", "correlation_yrna_pca" );
          }

          push( @report_files, "count_table_correlation",  "smallRNA_1mm_.+.snRNA.count.heatmap.png" );
          push( @report_files, "count_table_correlation",  "smallRNA_1mm_.+.snRNA.count.PCA.png" );
          push( @report_names, "correlation_snrna_heatmap", "correlation_snrna_pca" );

          push( @report_files, "count_table_correlation",  "smallRNA_1mm_.+.snoRNA.count.heatmap.png" );
          push( @report_files, "count_table_correlation",  "smallRNA_1mm_.+.snoRNA.count.PCA.png" );
          push( @report_names, "correlation_snorna_heatmap", "correlation_snorna_pca" );

          if ($hasGroupHeatmap) {
            push( @report_files, "count_table_correlation",        "smallRNA_1mm_.+.tRNA.count.Group.heatmap.png" );
            push( @report_files, "count_table_correlation",        "smallRNA_1mm_.+.tRNA.count.Group.Correlation.Cluster.png" );
            push( @report_names, "correlation_trna_group_heatmap", "correlation_trna_corr_cluster" );
          }
        }
      }

      if ( defined $config->{bowtie1_bacteria_group1_pm_table} ) {
        push( @report_files, "count_table_correlation",    "bacteria_group1_.*.Species.count.heatmap.png" );
        push( @report_files, "count_table_correlation",    "bacteria_group1_.*.Species.count.PCA.png" );
        push( @report_names, "correlation_group1_heatmap", "correlation_group1_pca" );

        if ($hasGroupHeatmap) {
          push( @report_files, "count_table_correlation",          "bacteria_group1_.*.Species.count.Group.heatmap.png" );
          push( @report_files, "count_table_correlation",          "bacteria_group1_.*.Species.count.Group.Correlation.Cluster.png" );
          push( @report_names, "correlation_group1_group_heatmap", "correlation_group1_corr_cluster" );
        }
      }

      if ( defined $config->{bowtie1_bacteria_group2_pm_table} ) {
        push( @report_files, "count_table_correlation",    "bacteria_group2_.*.Species.count.heatmap.png" );
        push( @report_files, "count_table_correlation",    "bacteria_group2_.*.Species.count.PCA.png" );
        push( @report_names, "correlation_group2_heatmap", "correlation_group2_pca" );

        if ($hasGroupHeatmap) {
          push( @report_files, "count_table_correlation",          "bacteria_group2_.*.Species.count.Group.heatmap.png" );
          push( @report_files, "count_table_correlation",          "bacteria_group2_.*.Species.count.Group.Correlation.Cluster.png" );
          push( @report_names, "correlation_group2_group_heatmap", "correlation_group2_corr_cluster" );
        }
      }

      if ( defined $config->{bowtie1_virus_group6_pm_table} ) {
        push( @report_files, "count_table_correlation",    "virus_group6_.*.Species.count.heatmap.png" );
        push( @report_files, "count_table_correlation",    "virus_group6_.*.Species.count.PCA.png" );
        push( @report_names, "correlation_group6_heatmap", "correlation_group6_pca" );

        if ($hasGroupHeatmap) {
          push( @report_files, "count_table_correlation",          "virus_group6_.*.Species.count.Group.heatmap.png" );
          push( @report_files, "count_table_correlation",          "virus_group6_.*.Species.count.Group.Correlation.Cluster.png" );
          push( @report_names, "correlation_group6_group_heatmap", "correlation_group6_corr_cluster" );
        }
      }

      if ( defined $config->{bowtie1_tRNA_pm_table} ) {
        push( @report_files, "count_table_correlation",     "^.*tRNA_pm_${task_name}.count.heatmap.png" );
        push( @report_files, "count_table_correlation",     "^.*tRNA_pm_${task_name}.count.PCA.png" );
        push( @report_names, "correlation_trnalib_heatmap", "correlation_trnalib_pca" );

        if ($hasGroupHeatmap) {
          push( @report_files, "count_table_correlation",           "^.*tRNA_pm_${task_name}.count.Group.heatmap.png" );
          push( @report_files, "count_table_correlation",           "^.*tRNA_pm_${task_name}.count.Group.Correlation.Cluster.png" );
          push( @report_names, "correlation_trnalib_group_heatmap", "correlation_trnalib_corr_cluster" );
        }
      }

      if ( defined $config->{bowtie1_rRNA_pm_table} ) {
        push( @report_files, "count_table_correlation",     "rRNA_pm_${task_name}.count.heatmap.png" );
        push( @report_files, "count_table_correlation",     "rRNA_pm_${task_name}.count.PCA.png" );
        push( @report_names, "correlation_rrnalib_heatmap", "correlation_rrnalib_pca" );

        if ($hasGroupHeatmap) {
          push( @report_files, "count_table_correlation",           "rRNA_pm_${task_name}.count.Group.heatmap.png" );
          push( @report_files, "count_table_correlation",           "rRNA_pm_${task_name}.count.Group.Correlation.Cluster.png" );
          push( @report_names, "correlation_rrnalib_group_heatmap", "correlation_rrnalib_corr_cluster" );
        }
      }

      if ( defined $config->{refseq_bacteria_table} ) {
        my $levels = [ 'species', 'genus', 'family', 'order', 'class', 'phylum'];

        for my $level (@$levels){
          push( @report_files, "count_table_correlation",     "${task_name}.${level}.aggregated.count.heatmap.png" );
          push( @report_files, "count_table_correlation",     "${task_name}.${level}.aggregated.count.PCA.png" );
          push( @report_names, "correlation_refseq_bacteria_${level}_agg_heatmap", "correlation_refseq_bacteria_${level}_agg_pca" );

          if ($hasGroupHeatmap) {
            push( @report_files, "count_table_correlation",           "${task_name}.${level}.aggregated.count.Group.heatmap.png" );
            push( @report_files, "count_table_correlation",           "${task_name}.${level}.aggregated.count.Group.Correlation.Cluster.png" );
            push( @report_names, "correlation_refseq_bacteria_${level}_agg_group_heatmap", "correlation_refseq_bacteria_${level}_agg_corr_cluster" );
          }

          push( @report_files, "count_table_correlation",     "${task_name}.${level}.estimated.count.heatmap.png" );
          push( @report_files, "count_table_correlation",     "${task_name}.${level}.estimated.count.PCA.png" );
          push( @report_names, "correlation_refseq_bacteria_${level}_est_heatmap", "correlation_refseq_bacteria_${level}_est_pca" );

          if ($hasGroupHeatmap) {
            push( @report_files, "count_table_correlation",           "${task_name}.${level}.estimated.count.Group.heatmap.png" );
            push( @report_files, "count_table_correlation",           "${task_name}.${level}.estimated.count.Group.Correlation.Cluster.png" );
            push( @report_names, "correlation_refseq_bacteria_${level}_est_group_heatmap", "correlation_refseq_bacteria_${level}_est_corr_cluster" );
          }
        }
      }
    }
    
    if ( defined $config->{host_length_dist_category} ) {
      push( @report_files, "host_length_dist_category", ".png" );
      push( @report_names, "host_length_dist_category" );
    }

    if ( defined $config->{nonhost_overlap_vis} ) {
      push( @report_files, "nonhost_overlap_vis", ".NonHost.Reads.Barplot.png" );
      push( @report_names, "nonhost_overlap_bar" );
    }

    if ( defined $config->{pairs} ) {
      if ( defined $config->{deseq2_host_genome_TotalReads_vis} ) {
        push( @report_files, "deseq2_host_genome_TotalReads_vis", ".DESeq2.Matrix.png" );
        push( @report_names, "deseq2_host_vis" );
      }elsif ( defined $config->{deseq2_host_genome_FeatureReads_vis} ) {
        push( @report_files, "deseq2_host_genome_FeatureReads_vis", ".DESeq2.Matrix.png" );
        push( @report_names, "deseq2_host_vis" );
      }elsif ( defined $config->{deseq2_host_genome_vis} ) {
        push( @report_files, "deseq2_host_genome_vis", ".DESeq2.Matrix.png" );
        push( @report_names, "deseq2_host_vis" );
      }
    }
    
    if(defined $config->{bacteria_count}){
      push( @copy_files, "bacteria_count", ".tsv\$", "bacteria_count", ".summary.png", "bacteria_count", ".summary.rpm.csv" );
      push( @report_files, "bacteria_count", ".summary.png",  "bacteria_count", ".summary.rpm.csv" );
      push( @report_names, "bacteria_count_vis", "bacteria_count_rpm" );
    }

    my $version_files = get_version_files($config);

    my $options = {
      "DE_fold_change" => [ getValue( $def, "DE_fold_change" ) ],
      "DE_pvalue"      => [ getValue( $def, "DE_pvalue" ) ],
      "DE_use_raw_pvalue"=> [ getValue( $def, "DE_use_raw_pvalue" ) ],
      "search_nonhost_genome" => [ getValue($def, "search_nonhost_genome") ],
      "normalize_by" => [ getValue($def, "normalize_by") ],
    };
    $config->{report} = {
      class                      => "CQS::BuildReport",
      perform                    => 1,
      target_dir                 => $def->{target_dir} . "/report",
      report_rmd_file            => "../Pipeline/SmallRNA.Rmd",
      additional_rmd_files       => "Functions.Rmd;../Pipeline/Pipeline.Rmd",
      parameterSampleFile1_ref   => \@report_files,
      parameterSampleFile1_names => \@report_names,
      parameterSampleFile2       => $options,
      parameterSampleFile3_ref   => \@copy_files,
      parameterSampleFile4       => $version_files,
      parameterSampleFile5       => $def->{software_version},
      parameterSampleFile6       => $def->{groups},
      sh_direct                  => 1,
      pbs                        => {
        "nodes"    => "1:ppn=1",
        "walltime" => "1",
        "mem"      => "10gb"
      },
    };
    push( @$summary_ref, "report" );
  }

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => $def->{target_dir} . "/sequencetask",
    option     => "",
    source     => {
      tasks => $individual_ref,
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
