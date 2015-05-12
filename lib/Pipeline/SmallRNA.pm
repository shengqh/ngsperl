#!/usr/bin/perl
package Pipeline::SmallRNA;

use strict;
use warnings;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performSmallRNA performSmallRNATask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub getSmallRNAConfig {
  my ($def) = @_;

  create_directory_or_die( $def->{target_dir} );

  my $cluster = $def->{cluster};
  if ( !defined $cluster ) {
    $cluster = "slurm";
  }

  my $fastq_remove_N = $def->{fastq_remove_N};
  my $run_cutadapt   = $def->{run_cutadapt};

  my $config = {
    general => {
      task_name => $def->{task_name},
      cluster   => $cluster
    },
    files => $def->{files}
  };

  my @individual = ();
  my @summary    = ();

  my $source_ref = "files";
  my $len_ref = "files";
  if ( !defined $fastq_remove_N || $fastq_remove_N ) {
    $config->{fastq_remove_N} = {
      class      => "CQS::FastqTrimmer",
      perform    => $fastq_remove_N,
      target_dir => $def->{target_dir} . "/fastq_remove_N",
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
    $len_ref = "fastq_remove_N";
    push @individual, "fastq_remove_N";
  }

  my $qc = {};
  if ( !defined $run_cutadapt || $run_cutadapt ) {
    my $adapter = $def->{adapter};
    if ( !defined $adapter ) {
      $adapter = "TGGAATTCTCGGGTGCCAAGG";
    }

    $qc = {
      fastqc_pre_trim => {
        class      => "QC::FastQC",
        perform    => 1,
        target_dir => $def->{target_dir} . "/fastqc_pre_trim",
        option     => "",
        source_ref => $source_ref,
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      },
      fastqc_pre_trim_summary => {
        class      => "QC::FastQCSummary",
        perform    => 1,
        target_dir => $def->{target_dir} . "/fastqc_pre_trim",
        cqstools   => $def->{cqstools},
        option     => "",
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      },
      cutadapt => {
        class      => "Cutadapt",
        perform    => 1,
        target_dir => $def->{target_dir} . "/cutadapt",
        option     => "-O 10 -m " . $def->{min_read_length},
        source_ref => $source_ref,
        adapter    => $adapter,
        extension  => "_clipped.fastq",
        sh_direct  => 1,
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "20gb"
        },
      },
      fastqc_post_trim => {
        class      => "QC::FastQC",
        perform    => 1,
        target_dir => $def->{target_dir} . "/fastqc_post_trim",
        option     => "",
        source_ref => [ "cutadapt", ".fastq.gz" ],
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      },
      fastqc_post_trim_summary => {
        class      => "QC::FastQCSummary",
        perform    => 1,
        target_dir => $def->{target_dir} . "/fastqc_post_trim",
        cqstools   => $def->{cqstools},
        option     => "",
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      }
    };
    $source_ref = [ "cutadapt", ".fastq.gz" ];
    $len_ref = "cutadapt";
    push @individual, ( "fastqc_pre_trim", "cutadapt", "fastqc_post_trim" );
    push @summary, ( "fastqc_pre_trim_summary", "fastqc_post_trim_summary" );
  }
  else {
    $qc = {
      fastqc => {
        class      => "QC::FastQC",
        perform    => 1,
        target_dir => $def->{target_dir} . "/fastqc",
        option     => "",
        source_ref => $source_ref,
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      },
      fastqc_summary => {
        class      => "QC::FastQCSummary",
        perform    => 1,
        target_dir => $def->{target_dir} . "/fastqc",
        cqstools   => $def->{cqstools},
        option     => "",
        cluster    => $cluster,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "2",
          "mem"      => "10gb"
        },
      }
    };
    push @individual, ("fastqc");
    push @summary,    ("fastqc_summary");
  }

  $config = merge( $config, $qc );

  #print Dumper($config);

  my $preparation = {
    fastq_len => {
      class      => "FastqLen",
      perform    => 1,
      target_dir => $def->{target_dir} . "/fastq_len",
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
    },
    identical => {
      class      => "FastqIdentical",
      perform    => 1,
      target_dir => $def->{target_dir} . "/identical",
      option     => "",
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
    identical_sequence_table => {
      class      => "CQS::SmallRNASequenceCountTable",
      perform    => 1,
      target_dir => $def->{target_dir} . "/identical_sequence_count_table",
      option     => "",
      source_ref => [ "identical", ".dupcount\$" ],
      cqs_tools  => $def->{cqstools},
      suffix     => "_sequence",
      sh_direct  => 1,
      cluster    => $cluster,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    identical_NTA => {
      class        => "CQS::FastqMirna",
      perform      => 1,
      target_dir   => $def->{target_dir} . "/identical_NTA",
      option       => "-l " . $def->{min_read_length},
      source_ref   => [ "identical", ".fastq.gz\$" ],
      seqcount_ref => [ "identical", ".dupcount\$" ],
      cqstools     => $def->{cqstools},
      extension    => "_clipped_identical_NTA.fastq.gz",
      sh_direct    => 1,
      cluster      => $cluster,
      pbs          => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },

    #1 mismatch search, NTA
    bowtie1_genome_1mm_NTA => {
      class         => "Bowtie1",
      perform       => 1,
      target_dir    => $def->{target_dir} . "/bowtie1_genome_1mm_NTA",
      option        => $def->{bowtie1_option_1mm},
      source_ref    => [ "identical_NTA", ".fastq.gz\$" ],
      bowtie1_index => $def->{bowtie1_index},
      samonly       => 0,
      sh_direct     => 1,
      cluster       => $cluster,
      pbs           => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=" . $def->{max_thread},
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },

    #not identical, for IGV
    bowtie1_genome_1mm_notidentical => {
      class         => "Bowtie1",
      perform       => 1,
      target_dir    => $def->{target_dir} . "/bowtie1_genome_1mm_notidentical",
      option        => $def->{bowtie1_option_1mm},
      source_ref    => $source_ref,
      bowtie1_index => $def->{bowtie1_index},
      samonly       => 0,
      sh_direct     => 0,
      cluster       => $cluster,
      pbs           => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=" . $def->{max_thread},
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
  };

  push @individual, ( "fastq_len", "identical", "identical_NTA", "bowtie1_genome_1mm_NTA", "bowtie1_genome_1mm_notidentical" );
  push @summary, ("identical_sequence_table");

  $config = merge( $config, $preparation );

  if ( defined $def->{coordinate} ) {
    push @individual, ( "bowtie1_genome_1mm_NTA_smallRNA_count", "bowtie1_genome_1mm_NTA_pmnames", "bowtie1_miRbase_pm", "bowtie1_miRbase_pm_count" );
    push @summary, ( "bowtie1_genome_1mm_NTA_smallRNA_table", "bowtie1_genome_1mm_NTA_smallRNA_category", "bowtie1_miRbase_pm_table" );

    my $count = {
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
        samtools        => $def->{samtools},
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

      #perfect match search to mirbase only
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
          "nodes"    => "1:ppn=" . $def->{max_thread},
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      bowtie1_miRbase_pm => {
        class         => "Alignment::Bowtie1",
        perform       => 1,
        target_dir    => $def->{target_dir} . "/bowtie1_miRbase_pm",
        option        => $def->{bowtie1_option_pm},
        source_ref    => [ "identical", ".fastq.gz\$" ],
        bowtie1_index => $def->{bowtie1_miRBase_index},
        samonly       => 0,
        sh_direct     => 1,
        cluster       => $cluster,
        pbs           => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=" . $def->{max_thread},
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      bowtie1_miRbase_pm_count => {
        class                   => "CQS::CQSChromosomeCount",
        perform                 => 1,
        target_dir              => $def->{target_dir} . "/bowtie1_miRbase_pm_count",
        option                  => $def->{mirbase_count_option},
        source_ref              => "bowtie1_miRbase_pm",
        seqcount_ref            => [ "identical", ".dupcount\$" ],
        perfect_mapped_name_ref => "bowtie1_genome_1mm_NTA_pmnames",
        cqs_tools               => $def->{cqstools},
        samtools                => $def->{samtools},
        sh_direct               => 1,
        cluster                 => $cluster,
        pbs                     => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "40gb"
        },
      },
      bowtie1_miRbase_pm_table => {
        class      => "CQS::CQSChromosomeTable",
        perform    => 1,
        target_dir => $def->{target_dir} . "/bowtie1_miRbase_pm_table",
        option     => "",
        source_ref => [ "bowtie1_miRbase_pm_count", ".xml" ],
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
      },
      sequencetask => {
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
      }
    };

    $config = merge( $config, $count );
  }

  return ($config);
}

sub performSmallRNA {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getSmallRNAConfig($def);

  if ($perform) {

    my $configFile = $def->{target_dir} . "/" . $def->{task_name} . ".config";
    open( SH, ">$configFile" ) or die "Cannot create $configFile";
    print SH Dumper($config);
    close(SH);

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
