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
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(initializeSmallRNADefaultOptions getSmallRNADefinition getPrepareConfig isVersion3 addNonhostDatabase addPositionVis addNonhostVis)] );

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
#  cqstools   => "/home/shengq1/cqstools/CQS.Tools.exe",
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

sub addNonhostDatabase {
  my ( $config, $def, $individual, $summary, $taskKey, $parentDir, $bowtieIndex, $sourceRef, $countOption, $tableOption, $count_ref ) = @_;

  my $bowtie1Task      = "bowtie1_" . $taskKey;
  my $bowtie1CountTask = "bowtie1_" . $taskKey . "_count";
  my $bowtie1TableTask = "bowtie1_" . $taskKey . "_table";

  if ( !defined $count_ref ) {
    $count_ref = [ "identical", ".dupcount\$" ];
  }

  addBowtie( $config, $def, $individual, $bowtie1Task, $parentDir, $bowtieIndex, $sourceRef, $def->{bowtie1_option_pm} );
  $config->{$bowtie1CountTask} = {
    class        => "CQS::CQSChromosomeCount",
    perform      => 1,
    target_dir   => $parentDir . "/" . $bowtie1CountTask,
    option       => $countOption,
    source_ref   => $bowtie1Task,
    seqcount_ref => $count_ref,
    cqs_tools    => $def->{cqstools},
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

  $config->{$bowtie1TableTask} = {
    class      => "CQS::CQSChromosomeTable",
    perform    => 1,
    target_dir => $parentDir . "/" . $bowtie1TableTask,
    option     => $tableOption,
    source_ref => [ $bowtie1CountTask, ".xml" ],
    cqs_tools  => $def->{cqstools},
    prefix     => $taskKey . "_",
    sh_direct  => 1,
    cluster    => $def->{cluster},
    pbs        => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };
  push @$individual, $bowtie1CountTask;
  push @$summary,    $bowtie1TableTask;
}

sub addPositionVis {
  my ( $config, $def, $summary, $taskName, $parentDir, $optionHash ) = @_;
  $config->{$taskName} = merge(
    $optionHash,
    {
      class                => "CQS::UniqueR",
      perform              => 1,
      target_dir           => $parentDir . "/$taskName",
      rtemplate            => "smallRnaPositionVis.R",
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

  $config->{$taskName} = merge(
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

sub initializeSmallRNADefaultOptions {
  my $def = shift;

  initDefaultValue( $def, "cluster",    "slurm" );
  initDefaultValue( $def, "max_thread", 8 );

  initDefaultValue( $def, "host_xml2bam",                    0 );
  initDefaultValue( $def, "bacteria_group1_count2bam",       0 );
  initDefaultValue( $def, "bacteria_group2_count2bam",       0 );
  initDefaultValue( $def, "fungus_group4_count2bam",         0 );
  initDefaultValue( $def, "host_bamplot",                    0 );
  initDefaultValue( $def, "read_correlation",                0 );
  initDefaultValue( $def, "perform_contig_analysis",         0 );
  initDefaultValue( $def, "perform_annotate_unmapped_reads", 0 );
  initDefaultValue( $def, "perform_nonhost_rRNA_coverage",   0 );
  initDefaultValue( $def, "perform_nonhost_tRNA_coverage",   0 );
  initDefaultValue( $def, "perform_host_rRNA_coverage",      0 );
  initDefaultValue( $def, "search_combined_nonhost",         0 );
  initDefaultValue( $def, "perform_report",                  0 );

  initDefaultValue( $def, "min_read_length",               16 );
  initDefaultValue( $def, "bowtie1_option_1mm",            "-a -m 100 --best --strata -v 1" );
  initDefaultValue( $def, "bowtie1_option_pm",             "-a -m 1000 --best --strata -v 0" );
  initDefaultValue( $def, "bowtie1_output_to_same_folder", 1 );
  initDefaultValue( $def, "fastq_remove_N",                1 );

  if ( defined $def->{run_cutadapt} && not defined $def->{perform_cutadapt} ) {
    $def->{perform_cutadapt} = $def->{run_cutadapt};
    $def->{run_cutadapt}     = undef;
  }
  initDefaultValue( $def, "perform_cutadapt", 1 );
  if ( $def->{perform_cutadapt} ) {
    initDefaultValue( $def, "adapter",         "TGGAATTCTCGGGTGCCAAGG" );
    initDefaultValue( $def, "cutadapt_option", "-m " . $def->{min_read_length} );
  }

  initDefaultValue( $def, "remove_sequences",          "" );
  initDefaultValue( $def, "fastq_remove_random",       0 );
  initDefaultValue( $def, "mirbase_count_option",      "-p hsa" );
  initDefaultValue( $def, "table_vis_group_text_size", 10 );
  initDefaultValue( $def, "sequencetask_run_time",     12 );

  initDefaultValue( $def, "DE_fold_change",              1.5 );
  initDefaultValue( $def, "DE_min_median_read_top",      2 );
  initDefaultValue( $def, "DE_min_median_read_smallRNA", 5 );
  initDefaultValue( $def, "DE_pvalue",                   0.05 );
  initDefaultValue( $def, "DE_use_raw_pvalue",           1 );
  initDefaultValue( $def, "DE_detected_in_both_group",   0 );
  initDefaultValue( $def, "DE_library_key",              "TotalReads" );

  initDefaultValue( $def, "smallrnacount_option",    "" );
  initDefaultValue( $def, "hasYRNA",                 0 );
  initDefaultValue( $def, "nonhost_table_option",    "--outputReadTable" );

  initDefaultValue( $def, "consider_miRNA_NTA", 1 );

  #search database
  initDefaultValue( $def, "search_not_identical",   1 );
  initDefaultValue( $def, "search_host_genome",     1 );
  initDefaultValue( $def, "search_nonhost_genome",  1 );
  initDefaultValue( $def, "search_nonhost_library", 1 );

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

  return $def;
}

sub getSmallRNADefinition {
  my ( $userdef, $genome ) = @_;

  my $def = merge( $userdef, $genome );

  $def = initializeSmallRNADefaultOptions($def);

  return $def;
}

sub getPrepareConfig {
  my ($def) = @_;

  #print Dumper($def);

  my $target_dir = create_directory_or_die( getValue( $def, "target_dir" ) );

  #check cqstools location
  my $cqstools = getValue( $def, "cqstools" );
  ( -e $cqstools ) or die "cqstools not exist: " . $cqstools;

  $def->{subdir} = 1;

  $def = initializeSmallRNADefaultOptions($def);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref ) = getPreprocessionConfig($def);

  my $class_independent_dir = create_directory_or_die( $target_dir . "/class_independent" );

  my $cluster = getValue( $def, "cluster" );

  #nta for microRNA and tRNA
  my $consider_miRNA_NTA = getValue( $def, "consider_miRNA_NTA" );
  my $consider_tRNA_NTA  = getValue( $def, "consider_tRNA_NTA" );

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
      cqstools   => $def->{cqstools},
      extension  => "_clipped_identical.fastq.gz",
      sh_direct  => 1,
      cluster    => $cluster,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
    identical_sequence_count_table => {
      class      => "CQS::SmallRNASequenceCountTable",
      perform    => 1,
      target_dir => $class_independent_dir . "/identical_sequence_count_table",
      option     => getValue( $def, "sequence_count_option" ),
      source_ref => [ "identical", ".dupcount\$" ],
      cqs_tools  => $def->{cqstools},
      suffix     => "_sequence",
      sh_direct  => 1,
      cluster    => $cluster,
      groups     => $def->{groups},
      pairs      => $def->{pairs},
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
  };

  push @$individual, ("identical");
  push @$summary,    ("identical_sequence_count_table");

  if ( $consider_miRNA_NTA && $consider_tRNA_NTA ) {
    $preparation->{identical_check_cca} = {
      class              => "SmallRNA::tRNACheckCCA",
      perform            => 1,
      target_dir         => $preprocessing_dir . "/identical_check_cca",
      option             => "",
      source_ref         => [ 'identical', '.fastq.gz$' ],
      untrimmedFastq_ref => $untrimed_ref,
      cqs_tools          => $def->{cqstools},
      sh_direct          => 1,
      pbs                => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb"
      },
    };
    push @$individual, ("identical_check_cca");
  }

  if ( $def->{special_sequence_file} ) {
    $config->{"special_sequence_count_table"} = {
      class                    => "CQS::ProgramWrapper",
      perform                  => 1,
      interpretor              => "python",
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

  if ( $consider_miRNA_NTA || $consider_tRNA_NTA ) {
    my $ccaaOption = $def->{consider_tRNA_NTA} ? "--ccaa" : "--no-ccaa";
    $preparation->{identical_NTA} = {
      class      => "SmallRNA::FastqSmallRnaNTA",
      perform    => 1,
      target_dir => $preprocessing_dir . "/identical_NTA",
      option     => $ccaaOption . " -l " . $def->{min_read_length},
      source_ref => [ "identical", ".fastq.gz\$" ],
      cqstools   => $def->{cqstools},
      extension  => "_clipped_identical_NTA.fastq.gz",
      sh_direct  => 1,
      cluster    => $cluster,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    };
    push @$individual, ("identical_NTA");
  }

  $config = merge( $config, $preparation );

  return ( $config, $individual, $summary, $cluster, $source_ref, $preprocessing_dir, $class_independent_dir );
}

1;
