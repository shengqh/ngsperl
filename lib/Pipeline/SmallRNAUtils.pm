#!/usr/bin/perl
package Pipeline::SmallRNAUtils;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(getSmallRNADefinition getPrepareConfig isVersion3 addNonhostDatabase addPositionVis addNonhostVis)] );

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
  my ( $config, $def, $individual, $summary, $taskKey, $parentDir, $bowtieIndex, $sourceRef, $countOption, $tableOption ) = @_;

  my $bowtie1Task      = "bowtie1_" . $taskKey;
  my $bowtie1CountTask = "bowtie1_" . $taskKey . "_count";
  my $bowtie1TableTask = "bowtie1_" . $taskKey . "_table";

  addBowtie( $config, $def, $individual, $bowtie1Task, $parentDir, $bowtieIndex, $sourceRef, $def->{bowtie1_option_pm} );
  $config->{$bowtie1CountTask} = {
    class        => "CQS::CQSChromosomeCount",
    perform      => 1,
    target_dir   => $parentDir . "/" . $bowtie1CountTask,
    option       => $countOption,
    source_ref   => $bowtie1Task,
    seqcount_ref => [ "identical", ".dupcount\$" ],
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
      sh_direct => 1,
      pbs       => {
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

sub initializeDefaultOptions {
  my $def = shift;

  initDefaultValue( $def, "cluster",            "slurm" );
  initDefaultValue( $def, "min_read_length",    16 );
  initDefaultValue( $def, "bowtie1_option_1mm", "-a -m 100 --best --strata -v 1" );
  initDefaultValue( $def, "bowtie1_option_pm",  "-a -m 1000 --best --strata -v 0" );
  initDefaultValue( $def, "fastq_remove_N",     1 );

  initDefaultValue( $def, "run_cutadapt", 1 );
  if ( $def->{run_cutadapt} ) {
    initDefaultValue( $def, "adapter",         "TGGAATTCTCGGGTGCCAAGG" );
    initDefaultValue( $def, "cutadapt_option", "-m " . $def->{min_read_length} );
  }

  initDefaultValue( $def, "remove_sequences",            "" );
  initDefaultValue( $def, "fastq_remove_random",         0 );
  initDefaultValue( $def, "mirbase_count_option",        "-p hsa" );
  initDefaultValue( $def, "table_vis_group_text_size",   10 );
  initDefaultValue( $def, "sequencetask_run_time",       12 );
  initDefaultValue( $def, "DE_show_gene_cluster",        1 );
  initDefaultValue( $def, "DE_pvalue",                   0.05 );
  initDefaultValue( $def, "DE_fold_change",              1.5 );
  initDefaultValue( $def, "DE_add_count_one",            0 );
  initDefaultValue( $def, "DE_min_median_read_top",      2 );
  initDefaultValue( $def, "DE_min_median_read_smallRNA", 5 );
  initDefaultValue( $def, "DE_top25only",                0 );
  initDefaultValue( $def, "DE_detected_in_both_group",   1 );
  initDefaultValue( $def, "DE_perform_wilcox",           0 );
  initDefaultValue( $def, "DE_use_raw_pvalue",           1 );
  initDefaultValue( $def, "DE_text_size",                10 );
  initDefaultValue( $def, "max_sequence_extension_base", 1 );
  initDefaultValue( $def, "top_read_number",             100 );

  initDefaultValue( $def, "perform_contig_analysis",     0 );
  initDefaultValue( $def, "smallrnacount_option",        "" );
  initDefaultValue( $def, "hasYRNA",                     0 );
  initDefaultValue( $def, "max_sequence_extension_base", "1" );
  initDefaultValue( $def, "nonhost_table_option",        "--outputReadTable" );

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

  my $additionalOption = $def->{hasYRNA} ? " --exportYRNA" : "";
  my $defaultOption = getValue( $def, "host_smallrnacount_option", "" );
  initDefaultValue( $def, "host_smallrnacount_option", $defaultOption . " --min_overlap 0.9 --offsets 0,1,2,-1,-2" . $additionalOption );

  $defaultOption = getValue( $def, "host_smallrnacounttable_option", "" );
  initDefaultValue( $def, "host_smallrnacounttable_option", $defaultOption . $additionalOption );

  return $def;
}

sub getSmallRNADefinition {
  my ( $userdef, $genome ) = @_;

  my $def = merge( $userdef, $genome );

  $def = initializeDefaultOptions($def);

  return $def;
}

sub getPrepareConfig {
  my ($def) = @_;

  #print Dumper($def);

  my $target_dir = create_directory_or_die( getValue( $def, "target_dir" ) );

  #check cqstools location
  my $cqstools = getValue( $def, "cqstools" );
  ( -e $cqstools ) or die "cqstools not exist: " . $cqstools;

  my $preprocessing_dir     = create_directory_or_die( $target_dir . "/preprocessing" );
  my $class_independent_dir = create_directory_or_die( $target_dir . "/class_independent" );

  $def = initializeDefaultOptions($def);
  my $cluster = getValue( $def, "cluster" );

  #remove contamination sequences from sequence kit before adapter trimming
  my $remove_sequences = getValue( $def, "remove_sequences" );

  #remove terminal N from fastq
  my $fastq_remove_N = getValue( $def, "fastq_remove_N" );

  #perform cutadapt to remove adapter
  my $run_cutadapt    = getValue( $def, "run_cutadapt" );
  my $cutadapt_option = getValue( $def, "cutadapt_option" );

  #for nextflex kit, we need to remove X bases after adapter trimming
  my $fastq_remove_random = getValue( $def, "fastq_remove_random" );

  #nta for microRNA and tRNA
  my $consider_miRNA_NTA = getValue( $def, "consider_miRNA_NTA" );
  my $consider_tRNA_NTA  = getValue( $def, "consider_tRNA_NTA" );

  my $max_sequence_extension_base = getValue( $def, "max_sequence_extension_base" );
  my $top_read_number             = getValue( $def, "top_read_number" );

  my $config = {
    general => {
      task_name => getValue( $def, "task_name" ),
      cluster   => $cluster
    },
    files => getValue( $def, "files" )
  };

  if ( defined $def->{groups} ) {
    $config->{groups} = $def->{groups};
  }

  if ( defined $def->{pairs} ) {
    $config->{pairs} = $def->{pairs};
  }

  my $individual = [];
  my $summary    = [];

  my $source_ref = "files";
  my $len_ref    = "files";
  if ( $fastq_remove_N && !$run_cutadapt ) {
    $config->{fastq_remove_N} = {
      class      => "CQS::FastqTrimmer",
      perform    => $fastq_remove_N,
      target_dir => $preprocessing_dir . "/fastq_remove_N",
      option     => "-n -z",
      extension  => "_trim.fastq.gz",
      source_ref => "files",
      cqstools   => $def->{cqstools},
      cluster    => $cluster,
      sh_direct  => 1,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "10gb"
      }
    };
    $source_ref = "fastq_remove_N";
    $len_ref    = "fastq_remove_N";
    push @$individual, "fastq_remove_N";
  }

  addFastQC( $config, $def, $individual, $summary, "fastqc_raw", $source_ref, $preprocessing_dir );

  if ( length($remove_sequences) ) {
    $config->{"remove_contamination_sequences"} = {
      class      => "CQS::Perl",
      perform    => 1,
      target_dir => $preprocessing_dir . "/remove_contamination_sequences",
      option     => $remove_sequences,
      output_ext => "_removeSeq.fastq.gz",
      perlFile   => "removeSequenceInFastq.pl",
      source_ref => $source_ref,
      sh_direct  => 0,
      cluster    => $cluster,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "20gb"
      },
    };
    push @$individual, ("remove_contamination_sequences");
    $source_ref = [ "remove_contamination_sequences", ".fastq.gz" ];
    $len_ref = "remove_contamination_sequences";

    addFastQC( $config, $def, $individual, $summary, "fastqc_post_remove", $source_ref, $preprocessing_dir );
  }

  if ($run_cutadapt) {
    my $adapter = $def->{adapter};
    if ( !defined $adapter ) {
      $adapter = "TGGAATTCTCGGGTGCCAAGG";
    }

    my $cutadapt_option = $def->{cutadapt_option};
    if ( !defined $cutadapt_option ) {
      $cutadapt_option = "-m " . $def->{min_read_length};
    }

    $config->{cutadapt} = {
      class                          => "Trimmer::Cutadapt",
      perform                        => 1,
      target_dir                     => $preprocessing_dir . "/cutadapt",
      option                         => $cutadapt_option,
      source_ref                     => $source_ref,
      adapter                        => $adapter,
      extension                      => "_clipped.fastq",
      random_bases_remove_after_trim => $fastq_remove_random,
      sh_direct                      => 0,
      cluster                        => $cluster,
      pbs                            => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    };
    push @$individual, "cutadapt";

    addFastQC( $config, $def, $individual, $summary, "fastqc_post_trim", [ "cutadapt", ".fastq.gz" ], $preprocessing_dir );
    $source_ref = [ "cutadapt", ".fastq.gz" ];
    $len_ref = "cutadapt";
  }

  my $fastqc_count_vis_files = undef;
  if ( length($remove_sequences) && $run_cutadapt ) {
    $fastqc_count_vis_files = {
      target_dir         => $preprocessing_dir . "/fastqc_post_trim",
      parameterFile2_ref => [ "fastqc_post_remove_summary", ".FastQC.summary.reads.tsv\$" ],
      parameterFile3_ref => [ "fastqc_post_trim_summary", ".FastQC.summary.reads.tsv\$" ],
    };
  }
  elsif ( length($remove_sequences) ) {
    $fastqc_count_vis_files = {
      target_dir         => $preprocessing_dir . "/fastqc_post_remove",
      parameterFile2_ref => [ "fastqc_post_remove_summary", ".FastQC.summary.reads.tsv\$" ],
    };
  }
  elsif ($run_cutadapt) {
    $fastqc_count_vis_files = {
      target_dir         => $preprocessing_dir . "/fastqc_post_trim",
      parameterFile2_ref => [ "fastqc_post_trim_summary", ".FastQC.summary.reads.tsv\$" ],
    };
  }

  if ( defined $fastqc_count_vis_files ) {
    $config->{"fastqc_count_vis"} = merge(
      {
        class              => "CQS::UniqueR",
        perform            => 1,
        rtemplate          => "countInFastQcVis.R",
        output_file        => ".countInFastQcVis.Result",
        output_file_ext    => ".Reads.csv",
        sh_direct          => 1,
        parameterFile1_ref => [ "fastqc_raw_summary", ".FastQC.summary.reads.tsv\$" ],
        pbs                => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },
      $fastqc_count_vis_files
    );
    push @$summary, ("fastqc_count_vis");
  }

  #print Dumper($config);
  $config->{"fastq_len"} = {
    class      => "CQS::FastqLen",
    perform    => 1,
    target_dir => $preprocessing_dir . "/fastq_len",
    option     => "",
    source_ref => $len_ref,
    cqstools   => $def->{cqstools},
    sh_direct  => 1,
    cluster    => $cluster,
    pbs        => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  };
  $config->{"fastq_len_vis"} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $preprocessing_dir . "/fastq_len",
    rtemplate                => "countTableVisFunctions.R,fastqLengthVis.R",
    output_file              => ".lengthDistribution",
    output_file_ext          => ".csv",
    parameterSampleFile1_ref => [ "fastq_len", ".len\$" ],
    parameterSampleFile2     => $def->{groups},
    sh_direct                => 1,
    pbs                      => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "10gb"
    },
  };
  push @$individual, ("fastq_len");
  push @$summary,    ("fastq_len_vis");

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
      option     => "--maxExtensionBase $max_sequence_extension_base -n $top_read_number --exportFastaNumber $top_read_number",
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
      untrimmedFastq_ref => "files",
      cqs_tools          => $def->{cqstools},
      sh_direct          => 0,
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
