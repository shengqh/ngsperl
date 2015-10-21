#!/usr/bin/perl
package Pipeline::ParclipSmallRNA;

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

our %EXPORT_TAGS = ( 'all' => [qw(getParclipSmallRNAConfig performParclipSmallRNA performParclipSmallRNATask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub getParclipSmallRNAConfig {
  my ($def) = @_;

  my ( $config, $individual_ref, $summary_ref, $cluster ) = getPrepareConfig( $def, 1 );
  my @individual = @{$individual_ref};
  my @summary    = @{$summary_ref};

  my $gsnap = {
    gsnap => {
      class                 => 'Alignment::Gsnap',
      perform               => 1,
      target_dir            => $def->{target_dir} . '/gsnap',
      option                => '-y 0 -z 0 -Y 0 -Z 0 -m 1 -Q --trim-mismatch-score 0 --trim-indel-score 0 --mode ttoc-nonstranded --gunzip',
      gsnap_index_directory => $def->{gsnap_index_directory},
      gsnap_index_name      => $def->{gsnap_index_name},
      source_ref            => [ 'identical_NTA', '.fastq.gz$' ],
      sh_direct             => 0,
      cluster               => $cluster,
      pbs                   => {
        'email'    => $def->{email},
        'nodes'    => '1:ppn=' . $def->{max_thread},
        'walltime' => '72',
        'mem'      => '40gb'
      }
    },
    gsnap_smallRNA_count => {
      class           => 'CQS::SmallRNACount',
      perform         => 1,
      target_dir      => $def->{target_dir} . '/gsnap_smallRNA_count',
      option          => '-s -e 4',
      source_ref      => 'gsnap',
      seqcount_ref    => [ 'identical_NTA', '.dupcount$' ],
      coordinate_file => $def->{coordinate},
      fasta_file      => $def->{coordinate_fasta},
      cqs_tools       => $def->{cqstools},
      sh_direct       => 0,
      cluster         => $cluster,
      pbs             => {
        'email'    => $def->{email},
        'walltime' => '72',
        'mem'      => '40gb',
        'nodes'    => '1:ppn=1'
      },
    },
    gsnap_smallRNA_t2c_summary => {
      class      => 'SmallRNA::T2CSummary',
      perform    => 1,
      target_dir => $def->{target_dir} . '/gsnap_smallRNA_t2c',
      option     => '',
      source_ref => [ 'gsnap_smallRNA_count', '.mapped.xml$' ],
      cqs_tools  => $def->{cqstools},
      sh_direct  => 0,
      cluster    => $cluster,
      pbs        => {
        'email'    => $def->{email},
        'walltime' => '72',
        'mem'      => '40gb',
        'nodes'    => '1:ppn=1'
      },
    },
    gsnap_smallRNA_table => {
      class      => "CQS::SmallRNATable",
      perform    => 1,
      target_dir => $def->{target_dir} . "/gsnap_smallRNA_table",
      option     => "",
      source_ref => [ "gsnap_smallRNA_count", ".mapped.xml" ],
      cqs_tools  => $def->{cqstools},
      prefix     => "smallRNA_parclip_",
      sh_direct  => 1,
      cluster    => $cluster,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "10",
        "mem"      => "10gb"
      },
    },
    gsnap_smallRNA_category => {
      class      => "CQS::SmallRNACategory",
      perform    => 1,
      target_dir => $def->{target_dir} . "/gsnap_smallRNA_category",
      option     => "",
      source_ref => [ "gsnap_smallRNA_count", ".info\$" ],
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

  push @individual, ( 'gsnap', 'gsnap_smallRNA_count' );
  push @summary, ( 'gsnap_smallRNA_t2c_summary', "gsnap_smallRNA_table", "gsnap_smallRNA_category" );

  $config = merge( $config, $gsnap );

  if ( defined $def->{search_3utr} && $def->{search_3utr} ) {
    ( defined $def->{utr3_db} ) or die "utr3_db should be defined with search_3utr for parclip data analysis.";
    ( -e $def->{utr3_db} ) or die "utr3_db defined but not exists : " . $def->{utr3_db};

    my $unmappedreads = {

      #extract unmapped reads
      unmappedReads => {
        class       => "CQS::Perl",
        perform     => 1,
        target_dir  => $def->{target_dir} . "/unmappedReads",
        perlFile    => "unmappedReadsToFastq.pl",
        source_ref  => [ "identical", ".fastq.gz\$" ],
        source2_ref => [ "gsnap_smallRNA_count", ".mapped.xml" ],
        output_ext  => "_clipped_identical.unmapped.fastq.gz",
        sh_direct   => 1,
        pbs         => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "10gb"
        },
      },

      #1 mismatch search
      unmappedReads_bowtie1_genome_1mm => {
        class         => "Bowtie1",
        perform       => 1,
        target_dir    => $def->{target_dir} . "/unmappedReads_bowtie1_genome_1mm",
        option        => $def->{bowtie1_option_1mm},
        source_ref    => [ "unmappedReads", ".fastq.gz\$" ],
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
      unmappedReads_bowtie1_genome_1mm_3utr_count => {
        class           => "CQS::SmallRNACount",
        perform         => 1,
        target_dir      => $def->{target_dir} . "/unmappedReads_bowtie1_genome_1mm_3utr_count",
        option          => "-m 0",
        source_ref      => [ "unmappedReads_bowtie1_genome_1mm", ".bam\$" ],
        fastq_files_ref => [ "unmappedReads", "fastq.gz\$" ],
        seqcount_ref    => [ "unmappedReads", ".dupcount\$" ],
        cqs_tools       => $def->{cqstools},
        coordinate_file => $def->{utr3_db},
        sh_direct       => 1,
        pbs             => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "20gb"
        },
      },

      #      unmappedReads_bowtie1_genome_1mm_3utr_count_target => {
      #        class        => "CQS::ParclipMirnaTarget",
      #        perform      => 1,
      #        target_dir   => $def->{target_dir} . "/unmappedReads_bowtie1_genome_1mm_3utr_count_target",
      #        option       => "",
      #        source_ref   => [ "t2c", ".xml\$" ],
      #        target_ref   => [ "utr3_count", ".xml\$" ],
      #        fasta_file   => $def->{fasta_file},
      #        refgene_file => $def->{refgene_file},
      #        cqs_tools    => $def->{cqstools},
      #        sh_direct    => 1,
      #        pbs          => {
      #          "email"    => $def->{email},
      #          "nodes"    => "1:ppn=1",
      #          "walltime" => "72",
      #          "mem"      => "20gb"
      #        },
      #      },
    };

    push( @individual, ( 'unmappedReads', 'unmappedReads_bowtie1_genome_1mm', 'unmappedReads_bowtie1_genome_1mm_3utr_count' ) );
    $config = merge( $config, $unmappedreads );
  }
  $config->{sequencetask} = {
    class      => 'CQS::SequenceTask',
    perform    => 1,
    target_dir => $def->{target_dir} . '/sequencetask',
    option     => '',
    source     => {
      step1 => \@individual,
      step2 => \@summary,
    },
    sh_direct => 0,
    cluster   => 'slurm',
    pbs       => {
      'email'    => $def->{email},
      'nodes'    => '1:ppn=' . $def->{max_thread},
      'walltime' => '72',
      'mem'      => '40gb'
    },
  };

  return ($config);
}

sub performParclipSmallRNA {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getParclipSmallRNAConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

sub performParclipSmallRNATask {
  my ( $def, $task ) = @_;

  my $config = getParclipSmallRNAConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
