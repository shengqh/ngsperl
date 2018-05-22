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

  my $groups            = $def->{groups};
  my $groups_vis_layout = $def->{groups_vis_layout};

  if ( defined $groups ) {
    if ( !defined $def->{tRNA_vis_group} ) {
      $def->{tRNA_vis_group} = $groups;
    }
  }

  my $t2c_dir = create_directory_or_die( $def->{target_dir} . "/t2c" );

  my $gsnap = {
    gsnap => {
      class                 => 'Alignment::Gsnap',
      perform               => 1,
      target_dir            => $t2c_dir . '/gsnap',
      option                => '-y 0 -z 0 -Y 0 -Z 0 -m 1 -Q --nofails --trim-mismatch-score 0 --trim-indel-score 0 --mode ttoc-nonstranded --gunzip',
      gsnap_index_directory => $def->{gsnap_index_directory},
      gsnap_index_name      => $def->{gsnap_index_name},
      source_ref            => [ 'identical_NTA', '.fastq.gz$' ],
      sh_direct             => 0,
      cluster               => $cluster,
      pbs                   => {
        'email'    => $def->{email},
        'nodes'    => '1:ppn=' . $def->{max_thread},
        'walltime' => '72',
        'mem'      => '80gb'
      }
    },
    gsnap_smallRNA_count => {
      class           => 'CQS::SmallRNACount',
      perform         => 1,
      target_dir      => $t2c_dir . '/gsnap_smallRNA_count',
      option          => '-s -e 4',
      source_ref      => 'gsnap',
      seqcount_ref    => [ 'identical', '.dupcount$' ],
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
    gsnap_smallRNA_table => {
      class      => "CQS::SmallRNATable",
      perform    => 1,
      target_dir => $t2c_dir . "/gsnap_smallRNA_table",
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
    gsnap_smallRNA_info => {
      class      => "CQS::CQSDatatable",
      perform    => 1,
      target_dir => $t2c_dir . "/gsnap_smallRNA_table",
      option     => "",
      source_ref => [ "gsnap_smallRNA_count", ".info" ],
      cqs_tools  => $def->{cqstools},
      prefix     => "smallRNA_parclip_",
      suffix     => ".mapped",
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

    gsnap_smallRNA_category => {
      class      => "CQS::SmallRNACategory",
      perform    => 1,
      target_dir => $t2c_dir . "/gsnap_smallRNA_category",
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
    gsnap_smallRNA_t2c => {
      class      => "CQS::ParclipT2CFinder",
      perform    => 1,
      target_dir => $t2c_dir . "/gsnap_smallRNA_t2c",
      option     => "-p 0.05 -e 0.013",
      source_ref => [ "gsnap_smallRNA_count", ".mapped.xml\$" ],
      cqs_tools  => $def->{cqstools},
      sh_direct  => 1,
      pbs        => {
        "email"    => $def->{email},
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "20gb"
      },
    },
    gsnap_smallRNA_t2c_summary => {
      class      => 'SmallRNA::T2CSummary',
      perform    => 1,
      target_dir => $t2c_dir . '/gsnap_smallRNA_t2c_table',
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
  };

  push @individual, ( 'gsnap', 'gsnap_smallRNA_count', 'gsnap_smallRNA_t2c' );
  push @summary, ( 'gsnap_smallRNA_table', 'gsnap_smallRNA_info', 'gsnap_smallRNA_category', 'gsnap_smallRNA_t2c_summary' );

  if ( defined $groups or defined $def->{tRNA_vis_group} ) {
    my $trna_vis_groups;
    if ( defined $def->{tRNA_vis_group} ) {
      $trna_vis_groups = $def->{tRNA_vis_group};
    }
    else {
      $trna_vis_groups = $groups;
    }

    addPositionVis(
      $config, $def,
      $summary_ref,
      "gsnap_tRNA_PositionVis",
      $t2c_dir,
      {
        output_file        => ".tRNAAnticodonPositionVis",
        output_file_ext    => ".tRNAAnticodonPositionVis.png",
        parameterFile1_ref => [ "gsnap_smallRNA_table", ".tRNA.count.position\$" ],
        parameterFile2_ref => [ "gsnap_smallRNA_info", ".mapped.count\$" ],
      }
    );
  }
  $config = merge( $config, $gsnap );

  if ( defined $def->{search_3utr} && $def->{search_3utr} ) {
    ( defined $def->{utr3_db} ) or die "utr3_db should be defined with search_3utr for parclip data analysis.";
    ( -e $def->{utr3_db} ) or die "utr3_db defined but not exists : " . $def->{utr3_db};
    my $utr3_dir = create_directory_or_die( $def->{target_dir} . "/3utr" );

    my $unmappedreads = {

      #extract unmapped reads
      unmappedReads => {
        class       => "CQS::Perl",
        perform     => 1,
        target_dir  => $utr3_dir . "/unmappedReads",
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
        class         => "Alignment::Bowtie1",
        perform       => 1,
        target_dir    => $utr3_dir . "/unmappedReads_bowtie1_genome_1mm",
        option        => $def->{bowtie1_option_1mm},
        source_ref    => [ "unmappedReads", ".fastq.gz\$" ],
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
      unmappedReads_bowtie1_genome_1mm_3utr_count => {
        class           => "CQS::SmallRNACount",
        perform         => 1,
        target_dir      => $utr3_dir . "/unmappedReads_bowtie1_genome_1mm_3utr_count",
        option          => "--noCategory",
        source_ref      => [ "unmappedReads_bowtie1_genome_1mm", ".bam\$" ],
        fastq_files_ref => [ "unmappedReads", ".fastq.gz\$" ],
        seqcount_ref    => [ "identical", ".dupcount\$" ],
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
      unmappedReads_bowtie1_genome_1mm_3utr_count_table => {
        class      => "CQS::CQSDatatable",
        perform    => 1,
        target_dir => $t2c_dir . "/unmappedReads_bowtie1_genome_1mm_3utr_count_table",
        option     => "-k 0 -i 1 -v 2 -e --fillMissingWithZero",
        source_ref => "unmappedReads_bowtie1_genome_1mm_3utr_count",
        cqs_tools  => $def->{cqstools},
        sh_direct  => 1,
        pbs        => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "10",
          "mem"      => "10gb"
        },
      },
      unmappedReads_bowtie1_genome_1mm_3utr_count_target_t2c => {
        class        => "CQS::ParclipTarget",
        perform      => 1,
        target_dir   => $utr3_dir . "/unmappedReads_bowtie1_genome_1mm_3utr_count_target_t2c",
        option       => "",
        source_ref   => [ "gsnap_smallRNA_t2c", ".xml\$" ],
        target_ref   => [ "unmappedReads_bowtie1_genome_1mm_3utr_count", ".xml\$" ],
        fasta_file   => $def->{fasta_file},
        refgene_file => $def->{refgene_file},
        cqs_tools    => $def->{cqstools},
        sh_direct    => 1,
        pbs          => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "20gb"
        },
      },

      unmappedReads_bowtie1_genome_1mm_3utr_count_target_all => {
        class        => "CQS::ParclipTarget",
        perform      => 1,
        target_dir   => $utr3_dir . "/unmappedReads_bowtie1_genome_1mm_3utr_count_target_all",
        option       => "",
        source_ref   => [ "gsnap_smallRNA_count", ".mapped.xml\$" ],
        target_ref   => [ "unmappedReads_bowtie1_genome_1mm_3utr_count", ".xml\$" ],
        fasta_file   => $def->{fasta_file},
        refgene_file => $def->{refgene_file},
        cqs_tools    => $def->{cqstools},
        sh_direct    => 1,
        pbs          => {
          "email"    => $def->{email},
          "nodes"    => "1:ppn=1",
          "walltime" => "72",
          "mem"      => "20gb"
        },
      },
    };

    push( @individual, ( 'unmappedReads', 'unmappedReads_bowtie1_genome_1mm', 'unmappedReads_bowtie1_genome_1mm_3utr_count', "unmappedReads_bowtie1_genome_1mm_3utr_count_target_t2c", "unmappedReads_bowtie1_genome_1mm_3utr_count_target_all" ) );
    push( @summary, 'unmappedReads_bowtie1_genome_1mm_3utr_count_table');
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
