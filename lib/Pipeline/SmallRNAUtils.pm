#!/usr/bin/perl
package Pipeline::SmallRNAUtils;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use Data::Dumper;

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(initializeSmallRNADefaultOptions 
  getSmallRNADefinition 
  getPrepareConfig 
  isVersion3 
  addNonhostDatabase 
  addPositionVis 
  addNonhostVis
  add_search_fasta
  )] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

#an example of parameter userdef
#my $userdef = {
#
#  #General options
#  task_name  => "parclip_NIH",
#  email      => "quanhu.sheng\@vanderbilt.edu",
#  target_dir => "/scratch/cqs/shengq1/vickers/20150925_parclip_3018-KCV-15/",
#  max_thread => 8,
#  cluster    => "slurm",
#  search_not_identical => 0,
#
#  #Default software parameter (don't change it except you really know it)
#  fastq_remove_N => 0,
#  adapter        => "TGGAATTCTCGGGTGCCAAGG",
#
#
#  #Data
#  files => {
#    "3018-KCV-15-15" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-15_ATGTCA_L006_R1_001.fastq.gz"],
#    "3018-KCV-15-36" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-36_CCAACA_L006_R1_001.fastq.gz"],
#    "3018-KCV-15-37" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-37_CGGAAT_L006_R1_001.fastq.gz"],
#    "3018-KCV-15-46" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-46_TCCCGA_L006_R1_001.fastq.gz"],
#    "3018-KCV-15-47" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-47_TCGAAG_L006_R1_001.fastq.gz"],
#  },
#};
#
#an example of paramter $genome
#my $genome = {
#  #genome database
#  mirbase_count_option  => "-p hsa",
#  coordinate            => "/scratch/cqs/shengq1/references/smallrna/hg19_miRBase20_ucsc-tRNA_ensembl75.bed",
#  coordinate_fasta      => "/scratch/cqs/shengq1/references/smallrna/hg19_miRBase20_ucsc-tRNA_ensembl75.bed.fa",
#  bowtie1_index         => "/scratch/cqs/shengq1/references/hg19_16569_MT/bowtie_index_1.1.2/hg19_16569_MT",
#  bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase20/bowtie_index_1.1.1/mature.dna",
#  gsnap_index_directory => "/scratch/cqs/shengq1/references/hg19_16569_MT/gsnap_index_k14_2015-06-23/",
#  gsnap_index_name      => "hg19_16569_MT",
#  star_index_directory => "/scratch/cqs/shengq1/references/hg19_16569_MT/STAR_index_v37.75_2.4.2a_sjdb49"
#};

sub isVersion3 {
  my $def = shift;
  my $result = defined $def->{version} && $def->{version} == 3;
  return $result;
}


sub initializeSmallRNADefaultOptions {
  my $def = shift;

  initDefaultValue( $def, "cluster",    "slurm" );
  initDefaultValue( $def, "max_thread", 8 );

  initDefaultValue( $def, "use_intermediate_dir", 1 );
  initDefaultValue( $def, "use_intermediate_dir_aside", 0 );
  
  initDefaultValue( $def, "host_xml2bam",                      0 );
  initDefaultValue( $def, "bacteria_group1_count2bam",         0 );
  initDefaultValue( $def, "bacteria_group2_count2bam",         0 );
  initDefaultValue( $def, "fungus_group4_count2bam",           0 );
  initDefaultValue( $def, "host_bamplot",                      0 );
  initDefaultValue( $def, "read_correlation",                  1 );
  initDefaultValue( $def, "perform_contig_analysis",           0 );
  initDefaultValue( $def, "perform_annotate_unmapped_reads",   0 );
  initDefaultValue( $def, "perform_nonhost_rRNA_coverage",     0 );
  initDefaultValue( $def, "perform_nonhost_tRNA_coverage",     0 );
  initDefaultValue( $def, "perform_host_tRNA_start_position",  1 );
  initDefaultValue( $def, "perform_host_rRNA_coverage",        1 );
  initDefaultValue( $def, "perform_host_length_dist_category", 1 );
  initDefaultValue( $def, "search_combined_nonhost",           0 );
  initDefaultValue( $def, "perform_bacteria_count",            1 );
  initDefaultValue( $def, "perform_report",                    1 );

  initDefaultValue( $def, "perform_short_reads_deseq2",        1 );
  initDefaultValue( $def, "perform_short_reads_source",        1 );

  initDefaultValue( $def, "perform_nonhost_genome_count",      1 );
  initDefaultValue( $def, "nonhost_genome_count_no_virus",     0 );

  initDefaultValue( $def, "min_read_length",               16 );
  initDefaultValue( $def, "bowtie1_option_2mm",            "-a --best --strata -v 2" );
  initDefaultValue( $def, "bowtie1_option_1mm",            "-a --best -m 100 --strata -v 1" );
  initDefaultValue( $def, "bowtie1_option_pm",             "-a --best -m 1000 --strata -v 0" );
  initDefaultValue( $def, "bowtie1_output_to_same_folder", 1 );
  initDefaultValue( $def, "fastq_remove_N",                1 );

  if ( defined $def->{run_cutadapt} && not defined $def->{perform_cutadapt} ) {
    $def->{perform_cutadapt} = $def->{run_cutadapt};
    $def->{run_cutadapt}     = undef;
  }
  initDefaultValue( $def, "perform_cutadapt", 1 );
  initDefaultValue( $def, "cutadapt_thread", 8 );

  initDefaultValue( $def, "fastq_len", 1 );

  initDefaultValue( $def, "use_first_read_after_trim", 1);
  
  if ( $def->{perform_cutadapt} ) {
    initDefaultValue( $def, "adapter",         "TGGAATTCTCGGGTGCCAAGG" );
    initDefaultValue( $def, "cutadapt_option", "-m " . $def->{min_read_length} );
    initDefaultValue( $def, "trim_poly_atgc",  1 );
    initDefaultValue( $def, "trim_base_quality_after_adapter_trim",  0 );
  }

  initDefaultValue( $def, "remove_sequences",          "'CCACGTTCCCGTGG;ACAGTCCGACGATC'" );
  initDefaultValue( $def, "fastq_remove_random",       0 );
  initDefaultValue( $def, "mirbase_count_option",      "-p hsa" );
  initDefaultValue( $def, "table_vis_group_text_size", 10 );
  initDefaultValue( $def, "sequencetask_run_time",     12 );

  initDefaultValue( $def, "top25cv_in_hca",     0 );

  #"TotalReads" or "None"
  initDefaultValue( $def, "normalize_by", "TotalReads");

  initDefaultValue( $def, "DE_fold_change",              1.5 );
  initDefaultValue( $def, "DE_min_median_read_top",      2 );
  initDefaultValue( $def, "DE_min_median_read_smallRNA", 5 );
  initDefaultValue( $def, "DE_pvalue",                   0.05 );
  initDefaultValue( $def, "DE_use_raw_pvalue",           1 );
  initDefaultValue( $def, "DE_detected_in_both_group",   0 );
  initDefaultValue( $def, "DE_library_key",              getValue($def, "normalize_by"));
  initDefaultValue( $def, "DE_cooksCutoff",              "FALSE" );
  initDefaultValue( $def, "use_pearson_in_hca",          0 );

  initDefaultValue( $def, "smallrnacount_option", "" );
  initDefaultValue( $def, "hasMicroRNAOnly",              0 );
  initDefaultValue( $def, "hasYRNA",              0 );
  initDefaultValue( $def, "nonhost_table_option", "--outputReadTable" );

  initDefaultValue( $def, "consider_miRNA_NTA", 1 );

  if ($def->{search_nonhost_genome_custom_group}){
    my $custom_group_name = getValue($def, "nonhost_genome_custom_group_name");
    if (defined $def->{customed_db}) {
      if (defined $def->{customed_db}{$custom_group_name}){
        $def->{bowtie1_custom_group_index} = $def->{customed_db}{$custom_group_name}{bowtie1_index};
        $def->{custom_group_species_map} = $def->{customed_db}{$custom_group_name}{species_map};
        die "File not exists: " . $def->{custom_group_species_map}  if !-e $def->{custom_group_species_map};
      }
    }
  }

  #search database
  initDefaultValue( $def, "search_nonhost_genome_custom_group_only",        0 );
  if (getValue( $def, "search_nonhost_genome_custom_group_only")){
    $def->{"search_not_identical"} = 0;
    $def->{"search_host_genome"} = 0;
    $def->{"search_nonhost_genome"} = 1;
    $def->{"search_nonhost_genome_custom_group"} = 1;
    $def->{"search_nonhost_library"} = 0;
    $def->{"perform_class_independent_analysis"} = 0;
    $def->{"blast_unmapped_reads"} = 0;
    $def->{"perform_report"} = 0;
    $def->{perform_host_genome_reads_deseq2} = 0;
  }else{
    initDefaultValue( $def, "search_not_identical",               0 );
    initDefaultValue( $def, "search_host_genome",                 1 );
    initDefaultValue( $def, "host_remove_all_mapped_reads",       0 );
    initDefaultValue( $def, "search_nonhost_genome",              1 );
    initDefaultValue( $def, "search_nonhost_library",             1 );
    initDefaultValue( $def, "perform_class_independent_analysis", 1 );
    initDefaultValue( $def, "perform_host_genome_reads_deseq2", 0 );
  }

  #blastn
  initDefaultValue( $def, "blast_top_reads",      0 );
  initDefaultValue( $def, "blast_unmapped_reads", 0 );
  initDefaultValue( $def, "blast_localdb",        "" );

  initDefaultValue( $def, "consider_tRNA_NTA", 1 );

  my $additionalOption = "";
  if ( $def->{hasSnRNA} ) {
    $additionalOption = $additionalOption . " --exportSnRNA";
  }
  if ( $def->{hasSnoRNA} ) {
    $additionalOption = $additionalOption . " --exportSnoRNA";
  }
  if ( $def->{hasYRNA} ) {
    $additionalOption = $additionalOption . " --exportYRNA";
  }
  my $defaultOption = getValue( $def, "host_smallrnacount_option", "" );
  initDefaultValue( $def, "host_smallrnacount_option", $defaultOption . " --min_overlap 0.9 --offsets 0,1,2,-1,-2" . $additionalOption );

  $defaultOption = getValue( $def, "host_smallrnacounttable_option", "" );
  initDefaultValue( $def, "host_smallrnacounttable_option", $defaultOption . $additionalOption );

  initDefaultValue( $def, "export_contig_details",       1 );
  initDefaultValue( $def, "max_sequence_extension_base", 1 );
  initDefaultValue( $def, "top_read_number",             100 );
  my $defaultSequenceCountOption = "--maxExtensionBase " . getValue( $def, "max_sequence_extension_base" ) .    #base
    " -n " . getValue( $def, "top_read_number" ) .                                                              #number of sequence
    " --exportFastaNumber " . getValue( $def, "top_read_number" ) .                                             #fasta
    ( getValue( $def, "export_contig_details" ) ? " --exportContigDetails" : "" );                              #contig_detail
  initDefaultValue( $def, "sequence_count_option", $defaultSequenceCountOption );

  #visualization
  initDefaultValue( $def, "use_least_groups", 0 );

  initDefaultValue( $def, "perform_nonhost_mappedToHost", 0 );

  initDefaultValue( $def, "perform_host_tRH_analysis", 0 );

  initDefaultValue( $def, "perform_host_tRnaFragmentHalves_analysis", 0 );

  initDefaultValue( $def, "correlation_rcode", "" );
  my $str = $def->{correlation_rcode};
  $str =~ s/^\s+|\s+$//g;
  if ( $str ne "" and $str !~ ";\$" ) {
    $str = $str . ";";
  }
  $def->{correlation_rcode} = $str;

  initDefaultValue( $def, "perform_permanova", 1 );

  return $def;
}

sub addNonhostDatabase {
  my ( $config, $def, $individual, $summary, $taskKey, $parentDir, $bowtieIndex, $sourceRef, $countOption, $tableOption, $count_ref, $nonhostXml ) = @_;

  my $bowtie1Task      = "bowtie1_" . $taskKey;
  my $bowtie1CountTask = "bowtie1_" . $taskKey . "_count";
  my $bowtie1TableTask = "bowtie1_" . $taskKey . "_table";

  if ( !defined $count_ref ) {
    $count_ref = [ "identical", ".dupcount\$" ];
  }

  my $intermediate_dir = getIntermidiateDir($parentDir, $def);
  #die($intermediate_dir);

  addBowtie( $config, $def, $individual, $bowtie1Task, $parentDir, $bowtieIndex, $sourceRef, $def->{bowtie1_option_pm} );
  $config->{$bowtie1CountTask} = {
    class        => "CQS::CQSChromosomeCount",
    perform      => 1,
    target_dir   => $intermediate_dir . "/" . $bowtie1CountTask,
    option       => $countOption,
    source_ref   => $bowtie1Task,
    seqcount_ref => $count_ref,
    sh_direct    => 1,
    cluster      => $def->{cluster},
    pbs          => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "72",
      "mem"       => "40gb"
    },
  };

  push @$nonhostXml, ( $bowtie1CountTask, ".xml" );

  $config->{$bowtie1TableTask} = {
    class      => "CQS::CQSChromosomeTable",
    perform    => 1,
    target_dir => $parentDir . "/" . $bowtie1TableTask,
    option     => $tableOption,
    source_ref => [ $bowtie1CountTask, ".xml" ],
    prefix     => $taskKey . "_",
    sh_direct  => 1,
    cluster    => $def->{cluster},
    pbs        => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push @$individual, $bowtie1CountTask;
  push @$summary,    $bowtie1TableTask;

  if ( $def->{perform_nonhost_mappedToHost_individual} ) {
    my $bowtie1readTask = "bowtie1_" . $taskKey . "_mappedreads";
    $config->{$bowtie1readTask} = {
      class                    => "CQS::ProgramIndividualWrapper",
      perform                  => 1,
      target_dir               => "${intermediate_dir}/$bowtie1readTask",
      option                   => "",
      interpretor              => "python3",
      program                  => "../SmallRNA/nonhostXmlToFastq.py",
      source_arg               => "-i",
      source_ref               => [ $bowtie1CountTask, ".count.mapped.xml" ],
      parameterSampleFile2_arg => "-f",
      parameterSampleFile2_ref => $sourceRef,
      output_to_same_folder    => 1,
      output_arg               => "-o",
      output_ext               => ".fastq.gz",
      sh_direct                => 1,
      pbs                      => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "10gb"
      },
    };
    push @$individual, $bowtie1readTask;

    my $bowtie1readMapTask = "bowtie1_" . $taskKey . "_mappedreads_host";
    addBowtie( $config, $def, $individual, $bowtie1readMapTask, $parentDir, $def->{bowtie1_index}, [$bowtie1readTask], $def->{bowtie1_option_2mm} );

    my $bowtie1readMapMismatchTask = "bowtie1_" . $taskKey . "_mappedreads_host_mismatch_table";
    $config->{$bowtie1readMapMismatchTask} = {
      class                    => "CQS::ProgramWrapper",
      perform                  => 1,
      target_dir               => "${parentDir}/$bowtie1readMapMismatchTask",
      option                   => "-m 2",
      interpretor              => "python3",
      program                  => "../SmallRNA/bamMismatchTable.py",
      parameterSampleFile1_arg => "-i",
      parameterSampleFile1_ref => [ $bowtie1readMapTask, ".bam\$" ],
      parameterSampleFile2_arg => "-f",
      parameterSampleFile2_ref => [ $bowtie1readTask, ".fastq.gz\$" ],
      parameterSampleFile3_arg => "-c",
      parameterSampleFile3_ref => $count_ref,
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
    push @$summary, $bowtie1readMapMismatchTask;
  }
}

sub addPositionVis {
  my ( $config, $def, $summary, $taskName, $parentDir, $optionHash ) = @_;
  $config->{$taskName} = merge_hash_left_precedent(
    $optionHash,
    {
      class                => "CQS::UniqueR",
      perform              => 1,
      target_dir           => $parentDir . "/$taskName",
      rtemplate            => "countTableVisFunctions.R,smallRnaPositionVis.R",
      output_file_ext      => ".allPositionBar.png",
      parameterFile2_ref   => [ "bowtie1_genome_1mm_NTA_smallRNA_info", ".mapped.count\$" ],
      parameterSampleFile1 => $def->{tRNA_vis_group},
      parameterSampleFile2 => $def->{groups_vis_layout},
      sh_direct            => 1,
      pbs                  => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    }
  );
  push @$summary, $taskName;
}

sub addNonhostVis {
  my ( $config, $def, $summary, $taskName, $parentDir, $optionHash ) = @_;

  $config->{$taskName} = merge_hash_left_precedent(
    $optionHash,
    {
      class      => "CQS::UniqueR",
      perform    => 1,
      target_dir => $parentDir . "/" . $taskName,

      #      parameterSampleFile1Order => $def->{groups_order},
      parameterSampleFile1 => $def->{groups},
      parameterSampleFile2 => $def->{groups_vis_layout},
      parameterFile3_ref   => [ "fastqc_count_vis", ".Reads.csv\$" ],
      rCode                => 'maxCategory=NA;textSize=9;groupTextSize=' . $def->{table_vis_group_text_size} . ';',
      sh_direct            => 1,
      pbs                  => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    }
  );

  push @$summary, $taskName;
}

sub getSmallRNADefinition {
  my ( $userdef, $genome ) = @_;

  my $def = merge_hash_left_precedent( $userdef, $genome );

  $def = initializeSmallRNADefaultOptions($def);

  return $def;
}

sub getPrepareConfig {
  my ($def) = @_;

  #print Dumper($def);

  my $target_dir = create_directory_or_die( getValue( $def, "target_dir" ) );

  $def->{subdir} = 1;

  $def = initializeSmallRNADefaultOptions($def);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $intermediate_dir = getIntermidiateDir($preprocessing_dir, $def);

  #nta for microRNA and tRNA
  my $hasMicroRNAOnly = getValue( $def, "hasMicroRNAOnly", 0 );
  my $consider_miRNA_NTA = getValue( $def, "consider_miRNA_NTA" );
  my $consider_tRNA_NTA  = getValue( $def, "consider_tRNA_NTA" ) && (!$hasMicroRNAOnly);

  if ( defined $def->{groups} ) {
    $config->{groups} = $def->{groups};
  }

  if ( defined $def->{pairs} ) {
    $config->{pairs} = $def->{pairs};
  }

  my $preparation = {
    identical => {
      class      => "CQS::FastqIdentical",
      perform    => 1,
      target_dir => $preprocessing_dir . "/identical",
      option     => "-l " . $def->{min_read_length},
      source_ref => $source_ref,
      extension  => "_clipped_identical.fastq.gz",
      sh_direct  => 1,
      use_first_read_after_trim => $def->{use_first_read_after_trim},
      cluster    => $cluster,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
  };
  push @$individual, ("identical");
  my $identical_ref = [ 'identical', '.fastq.gz$' ];
  my $host_identical_ref = [ 'identical', '.fastq.gz$' ];

  my $run_cutadapt     = getValue( $def, "perform_cutadapt");
  if($run_cutadapt){
    my $perform_identical_short = getValue( $def, "perform_identical_short", 0);
    if($perform_identical_short){
      $preparation->{identical_short} = {
        class      => "CQS::FastqIdentical",
        perform    => 1,
        target_dir => $preprocessing_dir . "/identical_short",
        option     => "-l 0",
        source_ref => [$source_ref->[0], ".fastq.short.gz"],
        extension  => "_clipped_identical_short.fastq.gz",
        sh_direct  => 1,
        use_first_read_after_trim => $def->{use_first_read_after_trim},
        cluster    => $cluster,
        pbs        => {
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "20gb"
        },
      };
      push @$individual, ("identical_short");
    }
  }

  if ( $consider_tRNA_NTA ) {
    $preparation->{identical_check_cca} = {
      class              => "SmallRNA::tRNACheckCCA",
      perform            => 1,
      target_dir         => $intermediate_dir . "/identical_check_cca",
      option             => "",
      source_ref         => [ 'identical', '.fastq.gz$' ],
      untrimmedFastq_ref => $untrimed_ref,
      sh_direct          => 1,
      pbs                => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    };
    push @$individual, ("identical_check_cca");
  }

  my $class_independent_dir;
  my $perform_class_independent_analysis = getValue( $def, "perform_class_independent_analysis", 1 );
  if ($perform_class_independent_analysis) {
    $class_independent_dir = create_directory_or_die( $target_dir . "/class_independent" );
    $preparation->{identical_sequence_count_table} = {
      class         => "CQS::SmallRNASequenceCountTable",
      perform       => 1,
      target_dir    => $class_independent_dir . "/identical_sequence_count_table",
      option        => getValue( $def, "sequence_count_option" ),
      source_ref    => [ "identical", ".dupcount\$" ],
      suffix        => "_sequence",
      sh_direct     => 1,
      cluster       => $cluster,
      groups        => $def->{groups},
      deseq2_groups => $def->{deseq2_groups},
      pairs         => $def->{pairs},
      pbs           => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => getValue($def, "identical_sequence_count_table_mem", "10gb")
      },
    };
    push @$summary, ("identical_sequence_count_table");

    if ( $def->{special_sequence_file} ) {
      $config->{"special_sequence_count_table"} = {
        class                    => "CQS::ProgramWrapper",
        perform                  => 1,
        interpretor              => "python3",
        program                  => "../SmallRNA/findSequence.py",
        target_dir               => $class_independent_dir . "/special_sequence_count_table",
        option                   => "",
        parameterSampleFile1_ref => [ "identical", ".dupcount\$" ],
        parameterFile1           => $def->{special_sequence_file},
        sh_direct                => 1,
        output_ext               => ".special_sequence.tsv",
        pbs                      => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "10",
          "mem"      => "10gb"
        },
      };
      push @$summary, ("special_sequence_count_table");
    }
  }

  if ( $consider_miRNA_NTA || $consider_tRNA_NTA ) {
    my $ccaaOption = $consider_tRNA_NTA ? "--ccaa" : "--no-ccaa";
    $preparation->{identical_NTA} = {
      class      => "SmallRNA::FastqSmallRnaNTA",
      perform    => 1,
      target_dir => $intermediate_dir . "/identical_NTA",
      option     => $ccaaOption . " -l " . $def->{min_read_length},
      source_ref => [ "identical", ".fastq.gz\$" ],
      extension  => "_clipped_identical_NTA.fastq.gz",
      sh_direct  => 1,
      cluster    => $cluster,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "20gb"
      },
    };
    push @$individual, ("identical_NTA");
    $host_identical_ref = [ "identical_NTA", ".fastq.gz\$" ];
  }

  $config = merge_hash_right_precedent( $config, $preparation );

  return ( $config, $individual, $summary, $cluster, $source_ref, $preprocessing_dir, $class_independent_dir, $identical_ref, $host_identical_ref );
}

sub add_chromosome_count {
  my ($config, $def, $tasks, $parent_dir, $count_task, $countOption, $bowtie1Task, $count_ref) = @_;

  $config->{$count_task} = {
    class        => "CQS::CQSChromosomeCount",
    perform      => 1,
    target_dir   => $parent_dir . "/" . $count_task,
    option       => $countOption,
    source_ref   => $bowtie1Task,
    seqcount_ref => $count_ref,
    sh_direct    => 1,
    pbs          => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push(@$tasks, $count_task);
}

sub add_absolute_chromosome_count {
  my ($config, $def, $tasks, $parent_dir, $count_task, $countOption, $bowtie1Task, $count_ref) = @_;

  $config->{$count_task} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform      => 1,
    target_dir   => $parent_dir . "/" . $count_task,
    option       => " -i __FILE__ -c __FILE2__ -o __NAME__.count $countOption

#__OUTPUT__

",
    source_ref   => [$bowtie1Task, '.bam$'],
    source_arg            => "-i",
    parameterSampleFile2_ref => $count_ref,
    parameterSampleFile2_arg => "-c",
    interpretor           => "python3",
    check_program         => 1,
    program               => "../SmallRNA/absolute_chromosome_count.py",
    output_to_same_folder => 1,
    output_arg            => "-o",
    output_to_folder      => 1,
    output_file_prefix    => "",
    output_file_ext       => ".count,.count.obj",
    sh_direct    => 1,
    pbs          => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "40gb"
    },
  };
  push(@$tasks, $count_task);
}

sub add_search_fasta {
  my ( $config, $def, $individual, $summary, $taskKey, $parentDir, $bowtie_index_task, $sourceRef, $bowtieOption, $countOption, $tableOption, $count_ref ) = @_;

  my $bowtie1Task      = "bowtie1_" . $taskKey . "_01_align";
  my $bowtie1CountTask = "bowtie1_" . $taskKey . "_02_count";
  my $bowtie1TableTask = "bowtie1_" . $taskKey . "_03_table";

  if ( !defined $count_ref ) {
    $count_ref = [ "identical", ".dupcount\$" ];
  }

  addBowtie( $config, $def, $individual, $bowtie1Task, $parentDir, $bowtie_index_task, $sourceRef, $bowtieOption );

  my $count_option = getValue($def, "search_fasta_count_option", "");

  add_absolute_chromosome_count($config, $def, $individual, $parentDir, $bowtie1CountTask, $count_option, $bowtie1Task, $count_ref);

  $config->{$bowtie1TableTask} = {
    class                    => "CQS::ProgramWrapper",
    perform                  => 1,
    interpretor              => "python3",
    program                  => "../SmallRNA/absolute_chromosome_count_table.py",
    target_dir               => $parentDir . "/$bowtie1TableTask",
    option                   => "",
    parameterSampleFile1_arg => "-i",
    parameterSampleFile1_ref => [$bowtie1CountTask, '.count$'],
    parameterSampleFile2_arg => "-s",
    parameterSampleFile2_ref => [$bowtie1CountTask, '.count.obj$'],
    sh_direct                => 1,
    output_arg               => "-o",
    output_ext               => ".count.txt",
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  };
  push (@$summary, $bowtie1TableTask);

  return($bowtie1Task, $bowtie1CountTask, $bowtie1TableTask);
}

1;
