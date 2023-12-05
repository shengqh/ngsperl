#!/usr/bin/perl
package Pipeline::Mixcr;

use strict;
use warnings;
use List::Util qw(first);
use File::Basename;
use Storable qw(dclone);
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

our %EXPORT_TAGS = ( 'all' => [qw(performMixcr)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );

  initDefaultValue( $def, "perform_preprocessing", 1 );

  initDefaultValue( $def, "perform_cutadapt", 0 );
  if ( $def->{perform_cutadapt} ) {
    #initDefaultValue( $def, "adapter", "CTGTCTCTTATA" );
    initDefaultValue( $def, "min_read_length", 30 );
    initDefaultValue( $def, "cutadapt_option", "-m " . $def->{min_read_length} );
    initDefaultValue( $def, "trim_polyA", 1 );
    initDefaultValue( $def, "trim_base_quality_after_adapter_trim",  0 );
  }

  initDefaultValue( $def, "max_thread",            8 );
  initDefaultValue( $def, "sequencetask_run_time", '24' );

  initDefaultValue( $def, "is_paired_end", 1 );

  initDefaultValue( $def, "perform_report", 1 );

  return $def;
}

sub getMixcrConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

  my $taskName = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $target_dir = $def->{target_dir};
  my $mixcr_species = getValue($def, "mixcr_species");
  my $mixcr_mem = getValue($def, "mixcr_mem", 40);
  my $mixcr_preset = getValue($def, "mixcr_preset");
  my $mixcr_thread = getValue($def, "mixcr_thread", 8);
  
  my $mixcr_analyze_rnaseq = "mixcr_analyze_rnaseq";
  $config->{$mixcr_analyze_rnaseq} = {
    class      => "CQS::ProgramWrapperOneToOne",
    perform    => 1,
    target_dir => "${target_dir}/$mixcr_analyze_rnaseq",
    option     => "

rm -f __NAME__.vdjca

mixcr analyze -Xmx${mixcr_mem}g \\
  --species $mixcr_species \\
  -t $mixcr_thread \\
  $mixcr_preset \\
  __FILE__ \\
  __NAME__

##__OUTPUT__

",
    program => "",
    check_program => 0,
    source_ref => $source_ref,
    source_arg => "",
    source_join_delimiter => " ",
    output_ext => ".clones_TRAD.tsv",
    sh_direct  => 0,
    cluster    => $cluster,
    output_to_same_folder => 0,
    pbs        => {
      "nodes"    => "1:ppn=1",
      "walltime" => "48",
      "mem"      => "${mixcr_mem}gb"
    },
  };

  if ( getValue( $def, "perform_report" ) ) {
    my @report_files = ();
    my @report_names = ();
    my @copy_files   = ();
    my $options      = {};

    my $version_files = get_version_files($config);

    if ( defined $config->{fastqc_raw_summary} ) {
      push( @report_files, "fastqc_raw_summary",                   ".FastQC.baseQuality.tsv.png" );
      push( @report_files, "fastqc_raw_summary",                   ".FastQC.sequenceGC.tsv.png" );
      push( @report_files, "fastqc_raw_summary",                   ".FastQC.adapter.tsv.png" );
      push( @report_names, "fastqc_raw_per_base_sequence_quality", "fastqc_raw_per_sequence_gc_content", "fastqc_raw_adapter_content" );
    }

    if ( defined $config->{fastqc_post_trim_summary} ) {
      push( @report_files, "fastqc_post_trim_summary",                   ".FastQC.baseQuality.tsv.png" );
      push( @report_files, "fastqc_post_trim_summary",                   ".FastQC.sequenceGC.tsv.png" );
      push( @report_files, "fastqc_post_trim_summary",                   ".FastQC.adapter.tsv.png" );
      push( @report_names, "fastqc_post_trim_per_base_sequence_quality", "fastqc_post_trim_per_sequence_gc_content", "fastqc_post_trim_adapter_content" );
    }

    if($def->{introduction_rmd}){
      $options->{introduction_rmd} = $def->{introduction_rmd};
    }

    $config->{report} = {
      class                      => "CQS::BuildReport",
      perform                    => 1,
      target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . "report",
      report_rmd_file            => "../Pipeline/RNASeq.Rmd",
      additional_rmd_files       => "../Pipeline/Pipeline.R;reportFunctions.R",
      parameterSampleFile1_ref   => \@report_files,
      parameterSampleFile1_names => \@report_names,
      parameterSampleFile2       => $options,
      parameterSampleFile3_ref   => \@copy_files,
      parameterSampleFile4       => $version_files,
      parameterSampleFile5       => $def->{software_version},
      parameterSampleFile6       => $def->{groups},
      
      sh_direct                  => 1,
      pbs                        => {
        "email"     => $def->{email},
        "emailType" => $def->{emailType},
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    };
    push( @$summary, "report" );
  }

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => $individual,
      step2 => $summary,
    },
    sh_direct => 0,
    cluster   => $cluster,
    pbs       => {
      "email"     => $def->{email},
      "emailType" => $def->{emailType},
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => $def->{sequencetask_run_time},
      "mem"       => "40gb"
    },
  };

  return ($config);
}

sub performMixcr {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getMixcrConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

1;
