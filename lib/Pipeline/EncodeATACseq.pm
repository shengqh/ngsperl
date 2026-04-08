#!/usr/bin/perl
package Pipeline::EncodeATACseq;

use strict;
use warnings;
use File::Basename;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::PeakPipelineUtils;
use Pipeline::Preprocession;
use Pipeline::WdlPipeline;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performEncodeATACseq performEncodeATACseqTask)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';


sub initializeDefaultOptions {
  my $def = shift;

  initDefaultValue( $def, "cluster",               "slurm" );
  initDefaultValue( $def, "max_thread",            8 );
  initDefaultValue( $def, "sequencetask_run_time", 12 );

  #Tasks
  initDefaultValue( $def, "sra_to_fastq", 0 );

  initDefaultValue( $def, "fastq_remove_N", 0 );

  initDefaultValue( $def, "perform_cutadapt", 0 );
  if ( getValue( $def, "perform_cutadapt" ) ) {
    #initDefaultValue( $def, "adapter", "CTGTCTCTTATA" );
    initDefaultValue( $def, "min_read_length", 30 );
    initDefaultValue( $def, "cutadapt_option", "-m " . $def->{min_read_length} );
    initDefaultValue( $def, "trim_polyA",      0 );
  } ## end if ( getValue( $def, "perform_cutadapt"...))
  else {
    #if no cutadapt, we need to call merge fastq to make sure the fastq files with the sample names user defined.
    $def->{merge_fastq} = 1;
  }

  if ( !defined $def->{"treatments"} ) {
    my $files  = getValue( $def, "files" );
    my $groups = {};
    for my $sample_name ( sort keys %$files ) {
      $groups->{$sample_name} = [$sample_name];
    }
    $def->{"treatments"} = $groups;
  } ## end if ( !defined $def->{"treatments"...})

  initDefaultValue( $def, "perform_croo_qc", 0 );

  initDefaultValue( $def, "perform_report", 0 );
  initDefaultValue( $def, "encode_option",  "" );
  initDefaultValue( $def, "caper_backend",  "local" );
  #$def->{encode_option} = "-b slurm"; #or empty to use local

  initDefaultValue( $def, "perform_diffbind", 0 );
  initDefaultValue( $def, "perform_homer",    0 );

  return $def;
} ## end sub initializeDefaultOptions


sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $tasks = [ @$individual, @$summary ];

  $def->{replicates} = $config->{groups};

  $def->{adapter} = "";

  my $croo_task     = addEncodeATACseq( $config, $def, $tasks, $target_dir, $source_ref, "encode_atacseq" );
  my $is_paired_end = is_paired_end($def);

  $config->{croo_bam} = {
    class         => "Encode::FindBamTask",
    source_ref    => [$croo_task],
    replicates    => $def->{groups},
    files_ref     => $untrimed_ref,
    has_cutadapt  => getValue( $def, "perform_cutadapt" ),
    is_paired_end => $is_paired_end,
  };

  $config->{croo_peak} = {
    class         => "Encode::FindPeakTask",
    source_ref    => [$croo_task],
    replicates    => $def->{groups},
    has_cutadapt  => getValue( $def, "perform_cutadapt" ),
    is_paired_end => $is_paired_end,
  };

  if ( $def->{perform_NFR_filter} ) {
    my $macs2_genome       = getValue( $def, "macs2_genome" );
    my $wdl                = $def->{"wdl"};
    my $encode_atac_folder = dirname( $wdl->{local}{encode_atacseq}{wdl_file} );
    print( $encode_atac_folder . "\n" );

    my $nfr_task = $croo_task . "_NFR_filter";
    $config->{$nfr_task} = {
      class         => "CQS::ProgramWrapperOneToOne",
      perform       => 1,
      target_dir    => $target_dir . "/" . $nfr_task,
      program       => "",
      check_program => 0,
      option        => "
echo samtools view -h -b -e 'tlen < 150 && tlen > -150' -o __NAME__.NFR_filtered.bam __FILE__ 
samtools view -h -b -e 'tlen < 150 && tlen > -150' -o __NAME__.NFR_filtered.bam __FILE__ 
samtools index __NAME__.NFR_filtered.bam

echo python3 $encode_atac_folder/src/encode_task_bam2ta.py __NAME__.NFR_filtered.bam --paired-end --mito-chr-name chrM --subsample 0 --mem-gb 6
python3 $encode_atac_folder/src/encode_task_bam2ta.py __NAME__.NFR_filtered.bam --paired-end --mito-chr-name chrM --subsample 0 --mem-gb 6

echo macs2 callpeak -t __NAME__.NFR_filtered.tn5.tagAlign.gz -f BED -n __NAME__.NFR_filtered -g $macs2_genome -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits
macs2 callpeak -t __NAME__.NFR_filtered.tn5.tagAlign.gz -f BED -n __NAME__.NFR_filtered -g $macs2_genome -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits

echo makeTagDirectory __NAME__.NFR_filtered.tagdir __NAME__.NFR_filtered.bam -format sam
makeTagDirectory __NAME__.NFR_filtered.tagdir __NAME__.NFR_filtered.bam -format sam

#__OUTPUT__ __FILE2__
",
      source_ref               => ["croo_bam"],
      output_file_ext          => ".NFR_filtered.bam,.NFR_filtered.tn5.tagAlign.gz",
      parameterSampleFile2_ref => $croo_task,                                          #make sure dependent on croo task.
      output_to_result_dir     => 1,
      output_to_same_folder    => 0,
      sh_direct                => 0,
      pbs                      => {
        "nodes"    => "1:ppn=8",
        "walltime" => "10",
        "mem"      => "10G"
      }
    };
    push( @$tasks, $nfr_task );
  } ## end if ( $def->{perform_NFR_filter...})

  if ( $def->{perform_diffbind} ) {
    my $bindName = $croo_task . "_diffbind";
    addDiffbind( $config, $def, $tasks, $target_dir, $bindName, "croo_bam", "croo_peak" );
    if ( getValue( $def, "perform_homer" ) ) {
      my $diffbind_homer = addHomerAnnotation( $config, $def, $tasks, $target_dir, $bindName, ".sig.bed" );
    }
  } ## end if ( $def->{perform_diffbind...})

  if ( $def->{perform_homer_analysis} ) {
    my $homer_tagdirectories_task = "homer_1_tagdirectories";
    add_homer_makeTagDirectory( $config, $def, $tasks, $target_dir, $homer_tagdirectories_task, "croo_bam", $croo_task );

    my $homer_mergePeaks_task = "homer_2_mergePeaks";
    add_homer_mergePeaks( $config, $def, $tasks, $target_dir, $homer_mergePeaks_task, "croo_peak", $croo_task );

    my $homer_annotatePeaks_task = "homer_3_annotatePeaks";
    add_homer_annotatePeaks( $config, $def, $tasks, $target_dir, $homer_annotatePeaks_task, $homer_tagdirectories_task, $homer_mergePeaks_task );

    if ( defined $def->{pairs} ) {
      my $deseq2taskname = addDEseq2( $config, $def, $tasks, "homer_peaks_raw_count", [ "homer_3_annotatePeaks", ".raw.count" ], $def->{target_dir}, $def->{DE_min_median_read} );
    }

  } ## end if ( $def->{perform_homer_analysis...})

  addSequenceTask( $config, $def, $tasks, $target_dir, $summary );

  return ($config);
} ## end sub getConfig


sub performEncodeATACseq {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getConfig($def);
  #print(Dumper($def));

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }
  return $config;
} ## end sub performEncodeATACseq


sub performEncodeATACseqTask {
  my ( $def, $task ) = @_;
  my $config = getConfig($def);

  performTask( $config, $task );

  return $config;
} ## end sub performEncodeATACseqTask

1;
