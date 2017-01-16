#!/usr/bin/perl
package Pipeline::ChIPSeq;

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

our %EXPORT_TAGS = ( 'all' => [qw(performChIPSeq performChIPSeqTask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;
  initValue( $def, "subdir", 0 );

  initValue( $def, "aligner", "bowtie1" );
  if ( $def->{aligner} eq "bowtie1" ) {
    initValue( $def, "bowtie1_option", "-v 1 -m 1 --best --strata" );
  }
  elsif ( $def->{aligner} eq "bowtie1" ) {
    initValue( $def, "bwa_option", "" );
  }

  initValue( $def, "peak_caller", "macs1" );
  if ( $def->{"peak_caller"} eq "macs1" ) {
    initValue( $def, "macs1_option", "-p 1e-9 -w -S --space=50" );
  }
  elsif ( $def->{peak_caller} eq "macs2" ) {
    my $macs2_genome = getValue( $def, "macs2_genome" );    #hs
    initValue( $def, "macs2_option", "--broad -B -q 0.01 -g " . $macs2_genome );
  }

  initValue( $def, "perform_rose",    0 );
  initValue( $def, "perform_bamplot", 0 );
  initValue( $def, "perform_chipqc",  0 );

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir ) = getPreprocessionConfig($def);

  my $email    = getValue( $def, "email" );
  my $cqstools = getValue( $def, "cqstools" );

  $config->{treatments} = getValue( $def, "treatments" );
  $config->{controls}   = getValue( $def, "controls" );

  if ( $def->{aligner} eq "bowtie1" ) {
    $config->{ $def->{aligner} } = {
      class         => "Alignment::Bowtie1",
      perform       => 1,
      target_dir    => "${target_dir}/" . $def->{aligner},
      option        => getValue( $def, "bowtie_option" ),
      fasta_file    => getValue( $def, "fasta_file" ),
      bowtie1_index => getValue( $def, "bowtie_index" ),
      source_ref    => $source_ref,
      sh_direct     => 0,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
  }
  elsif ( $def->{aligner} eq "bwa" ) {
    $config->{ $def->{aligner} } = {
      class      => "Alignment::BWA",
      perform    => 1,
      target_dir => "${target_dir}/" . $def->{aligner},
      option     => getValue( $def, "bwa_option" ),
      bwa_index  => getValue( $def, "bwa_index" ),
      source_ref => $source_ref,
      picard_jar => getValue( $def, "picard_jar" ),
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
  }
  else {
    die "Unknown alinger " . $def->{aligner};
  }

  push @$individual, ( $def->{aligner} );

  my $peakCallerTask = $def->{peak_caller} . "callpeak";
  if ( $def->{peak_caller} eq "macs1" ) {
    $config->{$peakCallerTask} = {
      class        => "Chipseq::MACS",
      perform      => 1,
      target_dir   => "${target_dir}/${peakCallerTask}",
      option       => getValue( $def, "macs1_option" ),
      source_ref   => [ $def->{aligner}, ".bam\$" ],
      groups_ref   => "treatments",
      controls_ref => "controls",
      sh_direct    => 0,
      pbs          => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
  }
  elsif ( $def->{peak_caller} eq "macs2" ) {
    $config->{$peakCallerTask} = {
      class        => "Chipseq::MACS2Callpeak",
      perform      => 1,
      target_dir   => "${target_dir}/$peakCallerTask",
      option       => getValue( $def, "macs2_option" ),
      source_ref   => [ $def->{aligner}, ".bam\$" ],
      groups_ref   => "treatments",
      controls_ref => "controls",
      sh_direct    => 0,
      pbs          => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
  }
  else {
    die "Unknown peak caller " . $def->{"peak_caller"};
  }
  push @$summary, ($peakCallerTask);

  if ( $def->{perform_rose} ) {
    my $roseTask = $def->{peak_caller} . "callpeak_bradner_rose2";
    $config->{$roseTask} = {
      class                => "Chipseq::BradnerRose2",
      perform              => 1,
      target_dir           => "${target_dir}/$roseTask",
      option               => "",
      source_ref           => $peakCallerTask,
      groups_ref           => "treatments",
      controls_ref         => "controls",
      pipeline_dir         => getValue( $def, "rose_folder" ),    #"/scratch/cqs/shengq1/local/bin/bradnerlab"
      binding_site_bed_ref => [ $peakCallerTask, ".bed\$" ],
      genome               => getValue( $def, "rose_genome" ),    #hg19,
      sh_direct            => 1,
      pbs                  => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    push @$summary, ($roseTask);
  }

  if ( getValue( $def, "perform_bamplot" ) ) {
    my $plotgroups = $def->{plotgroups};
    if ( !defined $plotgroups ) {
      my $files         = $def->{files};
      my @sortedSamples = sort keys %$files;
      $plotgroups = { getValue( $def, "task_name" ) => \@sortedSamples };
    }
    $config->{plotgroups} = $plotgroups;
    $config->{"bamplot"} = {
      class              => "Visualization::Bamplot",
      perform            => 1,
      target_dir         => "${target_dir}/bamplot",
      option             => getValue( $def, "bamplot_option" ),
      source_ref         => $def->{aligner},
      groups_ref         => "plotgroups",
      gff_file           => getValue( $def, "bamplot_gff" ),
      is_rainbow_color   => 0,
      is_draw_individual => 0,
      is_single_pdf      => 1,
      sh_direct          => 1,
      pbs                => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "1",
        "mem"      => "10gb"
      },
    };
    push @$summary, ("bamplot");
  }

  if ( getValue( $def, "perform_chipqc" ) ) {
    my $qctable = getValue( $def, "chipqc_table" );
    my $genome  = getValue( $def, "chipqc_genome" );    #hg19, check R ChIPQC package;
    $config->{chipqc} = {
      class      => "QC::ChipseqQC",
      perform    => 1,
      target_dir => "${target_dir}/chipqc",
      option     => "",
      source_ref => $def->{aligner},
      ,
      qctable       => $qctable,
      peaks_ref     => [ $peakCallerTask, ".bed\$" ],
      peak_software => "bed",
      genome        => $genome,
      sh_direct     => 0,
      pbs           => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    push @$summary, ("chipqc");
  }
  $config->{"sequencetask"} = {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_1 => $individual,
      step_2 => $summary,
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  return ($config);
}

sub performChIPSeq {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

sub performChIPSeqTask {
  my ( $def, $task ) = @_;
  my $config = getConfig($def);

  performTask( $config, $task );

  return $config;
}

1;
