#!/usr/bin/perl
package Pipeline::PipelineUtils;

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

our %EXPORT_TAGS =
  ( 'all' => [qw(getValue initPipelineOptions addPreprocess addFastQC addBlastn addBowtie addBamStat addDEseq2 addDeseq2Visualization addDeseq2SignificantSequenceBlastn getBatchGroups initDeseq2Options)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub getValue {
  my ( $def, $name, $defaultValue ) = @_;
  if ( defined $def->{$name} ) {
    return $def->{$name};
  }
  elsif ( defined $defaultValue ) {
    return $defaultValue;
  }
  else {
    die "Define $name in user definition first.";
  }
}

sub addFastQC {
  my ( $config, $def, $individual, $summary, $fastqcTask, $source_ref, $parentDir ) = @_;
  $config->{"$fastqcTask"} = {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => $parentDir . "/$fastqcTask",
    option     => "",
    source_ref => $source_ref,
    cluster    => $def->{cluster},
    pbs        => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  };

  my $summaryTask = $fastqcTask . "_summary";

  $config->{$summaryTask} = {
    class      => "QC::FastQCSummary",
    perform    => 1,
    target_dir => $parentDir . "/$fastqcTask",
    cqstools   => $def->{cqstools},
    option     => "",
    cluster    => $def->{cluster},
    ,
    pbs => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  };
  push @$individual, $fastqcTask;
  push @$summary,    $summaryTask;
}

sub addBlastn {
  my ( $config, $def, $summary, $blastTask, $fastaTask, $filePattern, $parentDir ) = @_;

  $config->{$blastTask} = {
    class      => "Blast::Blastn",
    perform    => 1,
    target_dir => $parentDir . "/$blastTask",
    option     => "",
    source_ref => [ $fastaTask, $filePattern ],
    sh_direct  => 0,
    localdb    => $def->{blast_localdb},
    cluster    => $def->{cluster},
    pbs        => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => "10",
      "mem"       => "10gb"
    }
  };

  push @$summary, $blastTask;
}

sub addBowtie {
  my ( $config, $def, $individual, $taskName, $parentDir, $bowtieIndex, $sourceRef, $bowtieOption ) = @_;

  $config->{$taskName} = {
    class         => "Alignment::Bowtie1",
    perform       => 1,
    target_dir    => $parentDir . "/" . $taskName,
    option        => $bowtieOption,
    source_ref    => $sourceRef,
    bowtie1_index => $bowtieIndex,
    samonly       => 0,
    sh_direct     => 1,
    mappedonly    => 1,
    cluster       => $def->{cluster},
    pbs           => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => "72",
      "mem"       => "40gb"
    },
  };

  push @$individual, $taskName;
}

sub addBamStat {
  my ( $config, $def, $summary, $taskName, $targetDir, $sourceRef ) = @_;

  $config->{$taskName} = {
    class                    => "CQS::UniqueR",
    target_dir               => $targetDir,
    perform                  => 1,
    rtemplate                => "../Samtools/BamStat.r",
    output_file              => ".bamstat.csv",
    sh_direct                => 1,
    parameterSampleFile1_ref => $sourceRef,
    cluster                  => $def->{cluster},
    pbs                      => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };
  push @$summary, $taskName;
}

sub addDEseq2 {
  my ( $config, $def, $summary, $taskKey, $countfileRef, $deseq2Dir, $DE_min_median_read ) = @_;

  my $taskName = "deseq2_" . $taskKey;

  $config->{$taskName} = {
    class                  => "Comparison::DESeq2",
    perform                => 1,
    target_dir             => $deseq2Dir . "/$taskName",
    option                 => "",
    source_ref             => "pairs",
    groups_ref             => "groups",
    countfile_ref          => $countfileRef,
    sh_direct              => 1,
    show_label_PCA         => $def->{show_label_PCA},
    use_pearson_in_hca     => $def->{use_pearson_in_hca},
    show_DE_gene_cluster   => $def->{DE_show_gene_cluster},
    pvalue                 => $def->{DE_pvalue},
    fold_change            => $def->{DE_fold_change},
    min_median_read        => $DE_min_median_read,
    add_count_one          => $def->{DE_add_count_one},
    top25only              => $def->{DE_top25only},
    detected_in_both_group => $def->{DE_detected_in_both_group},
    use_raw_p_value        => $def->{DE_use_raw_pvalue},
    text_size              => $def->{DE_text_size},
    cluster                => $def->{cluster},
    pbs                    => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "10",
      "mem"       => "10gb"
    },
  };
  push @$summary, $taskName;
  return $taskName;
}

sub addDeseq2Visualization {
  my ( $config, $def, $summary, $taskKey, $deseq2FileRef, $dataVisualizationDir, $layoutName ) = @_;

  my $taskName = "deseq2_" . $taskKey . "_vis";

  $config->{$taskName} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => $dataVisualizationDir . "/$taskName",
    rtemplate                => "DESeq2_all_vis.R",
    output_file              => ".${taskKey}.DESeq2.Matrix",
    output_file_ext          => ".png",
    parameterSampleFile1_ref => $deseq2FileRef,
    parameterSampleFile2     => $def->{$layoutName},
    rCode                    => 'useRawPvalue=' . $def->{DE_use_raw_pvalue} . ";",
    sh_direct                => 1,
    cluster                  => $def->{cluster},
    pbs                      => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "10gb"
    },
  };
  push @$summary, $taskName;
  return $taskName;
}

sub addDeseq2SignificantSequenceBlastn {
  my ( $config, $def, $summary, $deseq2Task, $parentDir ) = @_;

  my $fastaTask = $deseq2Task . "_sequences";
  $config->{$fastaTask} = {
    class      => "Blast::DESeq2SignificantReadToFasta",
    perform    => 1,
    target_dir => $parentDir . "/$fastaTask",
    option     => "",
    source_ref => [ $deseq2Task, "_DESeq2_sig.csv\$" ],
    sh_direct  => 1,
    cluster    => $def->{cluster},
    pbs        => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "10gb"
    }
  };

  push @$summary, ($fastaTask);

  addBlastn( $config, $def, $summary, $fastaTask . "_blastn", $fastaTask, ".fasta\$", $parentDir );
}

sub getBatchGroups {
  my ($def) = @_;
  my $files = $def->{files};
  my $result;

  if ( defined $def->{batch_groups_file_regex} ) {
    for my $regexName ( keys %{ $def->{batch_groups_file_regex} } ) {
      my $regex     = $def->{batch_groups_file_regex}->{$regexName};
      my $curGroups = {};
      $result->{$regexName} = $curGroups;
      for my $sample ( keys %$files ) {
        my $sampleFile = $files->{$sample}[0];
        my $group;
        if ( $sampleFile =~ /$regex/igs ) {
          my $groupName = $1;
          if ( !defined $curGroups->{$groupName} ) {
            $curGroups->{$groupName} = [];
          }
          my $groups = $curGroups->{$groupName};
          push @$groups, $sample;
        }
        else {
          die( $sample . " didn't match with regex " . $regex . " : " . $sampleFile );
        }
      }
    }
  }
  elsif ( defined $def->{batch_groups_name_regex} ) {
    for my $regexName ( keys %{ $def->{batch_groups_file_regex} } ) {
      my $regex     = $def->{batch_groups_name_regex}->{$regexName};
      my $curGroups = {};
      $result->{$regexName} = $curGroups;
      for my $sample ( keys %$files ) {
        my $group;
        if ( $sample =~ /$regex/igs ) {
          my $groupName = $1;
          if ( !defined $curGroups->{$groupName} ) {
            $curGroups->{$groupName} = [];
          }
          my $groups = $curGroups->{$groupName};
          push @$groups, $sample;
        }
        else {
          die( $sample . " didn't match with regex " . $regex );
        }
      }
    }
  }
  elsif ( defined $def->{batch_groups} ) {
    $result = $def->{batch_groups};
  }

  return $result;
}

sub initDeseq2Options {
  my $def = shift;
  initDefaultValue( $def, "DE_show_gene_cluster",      1 );
  initDefaultValue( $def, "DE_pvalue",                 0.05 );
  initDefaultValue( $def, "DE_fold_change",            2 );
  initDefaultValue( $def, "DE_add_count_one",          0 );
  initDefaultValue( $def, "DE_min_median_read",        5 );
  initDefaultValue( $def, "DE_top25only",              0 );
  initDefaultValue( $def, "DE_detected_in_both_group", 1 );
  initDefaultValue( $def, "DE_perform_wilcox",         0 );
  initDefaultValue( $def, "DE_use_raw_pvalue",         0 );
  initDefaultValue( $def, "DE_text_size",              10 );
  return $def;
}

1;
