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

our $VERSION = '0.01';

sub getSmallRNAConfig {
  my ($def) = @_;

  my ( $config, $individual_ref, $summary_ref, $cluster, $source_ref ) = getPrepareConfig($def, 1);
  my @individual = @{$individual_ref};
  my @summary    = @{$summary_ref};

  #print Dumper($config);

  my $bowtie1 = {

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

  push @individual, ( "bowtie1_genome_1mm_NTA", "bowtie1_genome_1mm_notidentical" );

  $config = merge( $config, $bowtie1 );

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
    saveConfig($def, $config);

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
