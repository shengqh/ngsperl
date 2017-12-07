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
  initDefaultValue( $def, "subdir", 0 );

  initDefaultValue( $def, "sra_to_fastq", 0 );

  my $pairend = getOption( $def, "pairend" );
  if ($pairend) {
    initDefaultValue( $def, "aligner", "bwa" );
  }
  else {
    initDefaultValue( $def, "aligner", "bowtie" );
  }
  if ( $def->{aligner} eq "bowtie1" ) {
    initDefaultValue( $def, "bowtie1_option", "-v 1 -m 1 --best --strata" );
  }
  elsif ( $def->{aligner} eq "bwa" ) {
    initDefaultValue( $def, "bwa_option", "" );
  }
  elsif ( $def->{aligner} eq "bowtie2" ) {
    initDefaultValue( $def, "bowtie2_option", "" );
  }

  initDefaultValue( $def, "peak_caller", "macs" );
  if ( $def->{"peak_caller"} eq "macs" ) {
    initDefaultValue( $def, "macs_option", "-p 1e-9 -w -S --space=50" );
  }
  elsif ( $def->{peak_caller} eq "macs2" ) {
    initDefaultValue( $def, "macs2_peak_type", "narrow" );
    if ( not defined $def->{"macs2_option"} ) {
      my $macs2_genome    = getValue( $def, "macs2_genome" );      #hzs
      my $macs2_peak_type = getValue( $def, "macs2_peak_type" );
      if ( $macs2_peak_type eq "narrow" ) {
        initDefaultValue( $def, "macs2_option", "-B -q 0.01 -g " . $macs2_genome );
      }
      else {
        initDefaultValue( $def, "macs2_option", "--broad -B -q 0.01 -g " . $macs2_genome );
      }
    }
  }

  initDefaultValue( $def, "perform_rose",     0 );
  initDefaultValue( $def, "perform_bamplot",  0 );
  initDefaultValue( $def, "perform_cleanbam", 0 );

  if ( $def->{perform_cleanbam} ) {
    initDefaultValue( $def, "minimum_maq",         10 );
    initDefaultValue( $def, "minimum_insert_size", 30 );
    initDefaultValue( $def, "maximum_insert_size", 1000 );
  }

  initDefaultValue( $def, "perform_chipqc",   0 );
  initDefaultValue( $def, "perform_diffbind", 0 );
  initDefaultValue( $def, "perform_enhancer", 0 );
  initDefaultValue( $def, "perform_multiqc",  1 );

  initDefaultValue( $def, "perform_enhancer", 0 );

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my $perform_chipqc = getValue( $def, "perform_chipqc" );

  my $perform_diffbind = getValue( $def, "perform_diffbind" );
  if ($perform_diffbind) {
    getValue( $def, "design_table" );
  }

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir ) = getPreprocessionConfig($def);
  my $step2 = [];

  my $email    = getValue( $def, "email" );
  my $cqstools = getValue( $def, "cqstools" );

  if ( $def->{aligner} eq "bowtie1" ) {
    $config->{bowtie1} = {
      class                 => "Alignment::Bowtie1",
      perform               => 1,
      target_dir            => "${target_dir}/" . getNextFolderIndex($def) . "bowtie1",
      option                => getValue( $def, "bowtie1_option" ),
      fasta_file            => getValue( $def, "bowtie1_fasta" ),
      bowtie1_index         => getValue( $def, "bowtie1_index" ),
      source_ref            => $source_ref,
      output_to_same_folder => 1,
      picard_jar            => getValue( $def, "picard_jar" ),
      mark_duplicates       => 1,
      sh_direct             => 0,
      pbs                   => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    $config->{bowtie1_summary} = {
      class                    => "CQS::UniqueR",
      perform                  => 1,
      rCode                    => "",
      target_dir               => $config->{bowtie1}->{target_dir},
      option                   => "",
      parameterSampleFile1_ref => [ "bowtie1", ".log" ],
      rtemplate                => "../Alignment/Bowtie1Summary.r",
      output_file              => $def->{task_name},
      output_file_ext          => ".csv",
      pbs                      => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb"
      },
    };

    push @$summary, ("bowtie1_summary");
  }
  elsif ( $def->{aligner} eq "bwa" ) {
    $config->{ $def->{aligner} } = {
      class                 => "Alignment::BWA",
      perform               => 1,
      target_dir            => "${target_dir}/" . getNextFolderIndex($def) . $def->{aligner},
      option                => getValue( $def, "bwa_option" ),
      bwa_index             => getValue( $def, "bwa_fasta" ),
      source_ref            => $source_ref,
      output_to_same_folder => 1,
      picard_jar            => getValue( $def, "picard_jar" ),
      mark_duplicates       => 1,
      sh_direct             => 0,
      pbs                   => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
  }
  elsif ( $def->{aligner} eq "bowtie2" ) {
    $config->{ $def->{aligner} } = {
      class                 => "Alignment::Bowtie2",
      perform               => 1,
      target_dir            => "${target_dir}/" . getNextFolderIndex($def) . $def->{aligner},
      option                => getValue( $def, "bowtie2_option" ),
      bowtie2_index         => getValue( $def, "bowtie2_index" ),
      source_ref            => $source_ref,
      output_to_same_folder => 1,
      picard_jar            => getValue( $def, "picard_jar" ),
      mark_duplicates       => 1,
      sh_direct             => 0,
      pbs                   => {
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
  my $bam_ref = [ $def->{aligner}, ".bam\$" ];

  push @$individual, ( $def->{aligner} );
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
      target_dir         => "${target_dir}/" . getNextFolderIndex($def) . "bamplot",
      option             => getValue( $def, "bamplot_option" ),
      source_ref         => $bam_ref,
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

  if ( $def->{perform_cleanbam} ) {
    my $taskName = $def->{aligner} . "_cleanbam";
    my $pairend = getValue( $def, "pairend" );
    addCleanBAM( $config, $def, $individual, $taskName, "${target_dir}/" . getNextFolderIndex($def) . $taskName, $bam_ref, $def->{pairend} );
    $bam_ref = [ $taskName, ".bam\$" ];
  }

  my $peakCallerTask;
  my $callFilePattern;

  if ( $def->{peak_caller} eq "macs" ) {
    $peakCallerTask            = "macs1callpeak";
    $callFilePattern           = ".name.bed\$";
    $config->{$peakCallerTask} = {
      class      => "Chipseq::MACS",
      perform    => 1,
      target_dir => "${target_dir}/" . getNextFolderIndex($def) . "${peakCallerTask}",
      option     => getValue( $def, "macs_option" ),
      source_ref => $bam_ref,
      groups     => $def->{"treatments"},
      controls   => $def->{"controls"},
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    $callFilePattern = ".name.bed\$";
  }
  elsif ( $def->{peak_caller} eq "macs2" ) {
    $peakCallerTask = "macs2callpeak";
    my $macs2option = getValue( $def, "macs2_option" );
    if ( $macs2option =~ /broad/ ) {
      $callFilePattern = "broadPeak.bed\$";
    }
    else {
      $callFilePattern = "narrowPeak.bed\$";
    }

    $peakCallerTask = $peakCallerTask . "_" . $def->{macs2_peak_type};
    $config->{$peakCallerTask} = {
      class      => "Chipseq::MACS2Callpeak",
      perform    => 1,
      target_dir => "${target_dir}/" . getNextFolderIndex($def) . "$peakCallerTask",
      option     => $macs2option,
      source_ref => $bam_ref,
      groups     => $def->{"treatments"},
      controls   => $def->{"controls"},
      sh_direct  => 0,
      pbs        => {
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
  push @$step2, ($peakCallerTask);
  if ( getValue( $def, "perform_homer_motifs" ) ) {
    addHomerMotif( $config, $def, $summary, $target_dir, $peakCallerTask, $callFilePattern );
  }

  if ( $def->{perform_rose} ) {
    my $roseTask = $peakCallerTask . "_bradner_rose2";
    $config->{$roseTask} = {
      class                => "Chipseq::Rose2",
      perform              => 1,
      target_dir           => "${target_dir}/" . getNextFolderIndex($def) . "$roseTask",
      option               => "",
      source_ref           => $bam_ref,
      groups               => $def->{"treatments"},
      controls             => $def->{"controls"},
      pipeline_dir         => getValue( $def, "rose_folder" ),                             #"/scratch/cqs/shengq1/local/bin/bradnerlab"
      binding_site_bed_ref => [ $peakCallerTask, ".bed\$" ],
      genome               => getValue( $def, "rose_genome" ),                             #hg19,
      sh_direct            => 1,
      pbs                  => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    push @$step2, ($roseTask);
  }

  if ($perform_chipqc) {
    my $genome = getValue( $def, "chipqc_genome" );      #hg19, check R ChIPQC package;
    my $chipqc_taskname = $peakCallerTask . "_chipqc";
    $config->{$chipqc_taskname} = {
      class          => "QC::ChipseqQC",
      perform        => 1,
      target_dir     => "${target_dir}/" . getNextFolderIndex($def) . $chipqc_taskname,
      option         => "",
      source_ref     => $bam_ref,
      groups         => $def->{"treatments"},
      controls       => $def->{"controls"},
      qctable        => $def->{"design_table"},
      peaks_ref      => [ $peakCallerTask, ".bed\$" ],
      peak_software  => "bed",
      genome         => $genome,
      combined       => getValue( $def, "chipqc_combined", 1 ),
      blacklist_file => $def->{"blacklist_file"},
      chromosomes    => $def->{"chipqc_chromosomes"},
      sh_direct      => 0,
      pbs            => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    push @$summary, ($chipqc_taskname);
  }

  if ($perform_diffbind) {
    my $bindName = $peakCallerTask . "_diffbind";
    $config->{$bindName} = {
      class                   => "Comparison::DiffBind",
      perform                 => 1,
      target_dir              => "${target_dir}/" . getNextFolderIndex($def) . "${bindName}",
      option                  => "",
      source_ref              => $bam_ref,
      groups                  => $def->{"treatments"},
      controls                => $def->{"controls"},
      design_table            => getValue( $def, "design_table" ),
      peaks_ref               => [ $peakCallerTask, ".bed\$" ],
      peak_software           => "bed",
      homer_annotation_genome => $def->{homer_annotation_genome},
      sh_direct               => 0,
      pbs                     => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    push @$summary, ($bindName);
  }

  if ( getValue( $def, "perform_enhancer" ) ) {
    addEnhancer( $config, $def, $individual, $summary, $target_dir, $peakCallerTask . "_enhancer", $bam_ref, [ $peakCallerTask, ".bed\$" ] );
  }

  if ( getValue( $def, "perform_multiqc" ) ) {
    addMultiQC( $config, $def, $summary, $target_dir, $target_dir );
  }
  $config->{"sequencetask"} = {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_1 => $individual,
      step_2 => $step2,
      step_3 => $summary,
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
