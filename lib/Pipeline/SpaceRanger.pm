#!/usr/bin/perl
package Pipeline::SpaceRanger;

use strict;
use warnings;
use List::Util qw(first);
use File::Basename;
use Storable qw(dclone);
use File::Slurp;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use Data::Dumper;
use Hash::Merge qw( merge );
use Storable    qw(dclone);
use scRNA::Modules;

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(
        add_copy_report
        add_SpaceRanger_summary
        performSpaceRanger)
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';


sub initializeSpaceRangerDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType",             "FAIL" );
  initDefaultValue( $def, "cluster",               "slurm" );
  initDefaultValue( $def, "perform_preprocessing", 0 );

  return $def;
} ## end sub initializeSpaceRangerDefaultOptions


sub add_copy_report {
  my ( $config, $def, $target_dir, $spaceranger_task, $tasks ) = @_;

  my $copy_report_task = "${spaceranger_task}_report";
  $config->{$copy_report_task} = {
    class         => "CQS::ProgramWrapperOneToOne",
    target_dir    => "$target_dir/$copy_report_task",
    perform       => 1,
    program       => "",
    check_program => 0,
    option        => "

cp __FILE__ __NAME__.web_summary.html

",
    parameterSampleFile1_ref => [ $spaceranger_task, "web_summary.html" ],
    sh_direct                => 1,
    no_docker                => 1,
    no_output                => 1,
    output_ext               => ".web_summary.html",
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "2gb"
    }
  };
  push( @$tasks, $copy_report_task );
  return ($copy_report_task);
} ## end sub add_copy_report


sub add_SpaceRanger_summary {
  my ( $config, $def, $target_dir, $spaceranger_task, $tasks, $copy_report_task ) = @_;

  my $spaceranger_summary_task = "${spaceranger_task}_summary";
  $config->{$spaceranger_summary_task} = {
    class                    => "CQS::UniqueRmd",
    target_dir               => "$target_dir/$spaceranger_summary_task",
    perform                  => 1,
    option                   => "",
    report_rmd_file          => "../scRNA/space_ranger_summary.Rmd",
    additional_rmd_files     => "../CQS/reportFunctions.R;../scRNA/scRNA_func.r;../scRNA/space_ranger_summary_cell.Rmd",
    parameterSampleFile1_ref => [ $spaceranger_task, "metrics_summary.csv" ],
    parameterSampleFile2     => {
      task_name   => getValue( $def, "task_name" ),
      email       => getValue( $def, "email" ),
      affiliation => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
    },
    parameterSampleFile3_ref => [ $copy_report_task, ".html" ],
    output_file_ext          => ".spaceranger.html",
    output_other_ext         => ".spaceranger.html",
    sh_direct                => 0,
    no_docker                => getValue( $def, "no_docker", 0 ),
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "8",
      "mem"      => "40gb"
    }
  };

  push( @$tasks, $spaceranger_summary_task );
  return $spaceranger_summary_task;
} ## end sub add_SpaceRanger_summary


sub getSpaceRangerConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeSpaceRangerDefaultOptions($def);

  my $project_name = $def->{task_name};

  my $email = $def->{email};

  # Disable preprocessing for spatial transcriptome
  $def->{perform_preprocessing} = 0;

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $tasks = [ @$individual, @$summary ];

  my $target_dir = getValue( $def, "target_dir" );

  if ( defined $def->{raw_spaceranger_metrics} ) {
    my $raw_spaceranger = {};
    for my $k ( keys %{ $def->{raw_spaceranger_metrics} } ) {
      my $metrics_file     = $def->{raw_spaceranger_metrics}->{$k}[0];
      my $web_summary_file = $metrics_file;
      $web_summary_file =~ s/metrics_summary\.csv/web_summary.html/;
      $raw_spaceranger->{$k} = [ $metrics_file, $web_summary_file ];
    } ## end for my $k ( keys %{ $def...})
    $config->{raw_spaceranger} = $raw_spaceranger;

    my $raw_spaceranger_copy_report_task = add_copy_report( $config, $def, $target_dir, "raw_spaceranger", $tasks );
    my $raw_spaceranger_summary_task     = add_SpaceRanger_summary( $config, $def, $target_dir, "raw_spaceranger", $tasks, $raw_spaceranger_copy_report_task );
  } ## end if ( defined $def->{raw_spaceranger_metrics...})

  my $spaceranger_task          = "spaceranger";
  my $fastq_folder              = getValue( $def, "fastq_folder" );
  my $spaceranger_transcriptome = getValue( $def, "spaceranger_transcriptome" );
  my $spaceranger_probe_set     = getValue( $def, "spaceranger_probe_set" );
  my $spaceranger_jobmode       = getValue( $def, "spaceranger_jobmode", "slurm" );

  $config->{files}                  = getValue( $def, "files" );
  $config->{H_E_brightfield_images} = getValue( $def, "H_E_brightfield_images" );
  $config->{manual_alignment_json}  = $def->{"manual_alignment_json"};

  my $json_option = "";
  if ( defined $def->{manual_alignment_json} ) {
    $json_option = "--loupe-alignment __FILE4__";
  }

  if ( defined $def->{id_to_samples} ) {
    $config->{id_to_samples} = getValue( $def, "id_to_samples" );
  }
  else {
    my $all_samples = [];
    for my $k ( keys %{ $config->{files} } ) {
      push( @$all_samples, $k );
    }
    $config->{id_to_samples} = { map { $_ => [$_] } @$all_samples };
  } ## end else [ if ( defined $def->{id_to_samples...})]

  $config->{$spaceranger_task} = {
    class         => "CQS::ProgramWrapperOneToOne",
    target_dir    => "$target_dir/$spaceranger_task",
    perform       => 1,
    program       => "",
    check_program => 0,
    option        => "

rm -rf __NAME__ ____NAME__.mro

spaceranger count --disable-ui \\
  --id __NAME__ \\
  --transcriptome $spaceranger_transcriptome \\
  --probe-set $spaceranger_probe_set \\
  --create-bam false \\
  --cytaimage __FILE__ \\
  --image __FILE2__ \\
  --sample __FILE3__ $json_option \\
  --fastqs $fastq_folder \\
  --jobmode $spaceranger_jobmode

rm -rf __NAME__/SPATIAL_RNA_COUNTER_CS \\
  __NAME__/extras \\
  __NAME__/_filelist \\
  __NAME__/_finalstate \\
  __NAME__/_invocation \\
  __NAME__/_jobmode \\
  __NAME__/_mrosource \\
  __NAME__/_perf \\
  __NAME__/_perf._truncated_ \\
  __NAME__/_sitecheck \\
  __NAME__/_tags \\
  __NAME__/_uuid \\
  __NAME__/_vdrkill 

",
    parameterSampleFile1_ref => "files",
    parameterSampleFile2_ref => "H_E_brightfield_images",
    parameterSampleFile3_ref => "id_to_samples",
    sh_direct                => 1,
    no_docker                => 1,
    no_output                => 1,
    output_ext               => "__NAME__/outs/web_summary.html,__NAME__/outs/metrics_summary.csv,__NAME__/outs/segmented_outputs/filtered_feature_cell_matrix.h5",
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "40gb"
    }
  };
  if ( defined $def->{manual_alignment_json} ) {
    $config->{$spaceranger_task}{parameterSampleFile4_ref} = "manual_alignment_json";
  }

  push( @$tasks, $spaceranger_task );

  my $spaceranger_copy_report_task = add_copy_report( $config, $def, $target_dir, $spaceranger_task, $tasks );
  my $spaceranger_summary_task     = add_SpaceRanger_summary( $config, $def, $target_dir, $spaceranger_task, $tasks, $spaceranger_copy_report_task );

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => { step1 => $tasks, },
    sh_direct  => 0,
    cluster    => $cluster,
    pbs        => {
      "nodes"    => "1:ppn=" . $def->{max_thread},
      "walltime" => $def->{sequencetask_run_time},
      "mem"      => "40gb"
    },
  };

  return ($config);
} ## end sub getSpaceRangerConfig


sub performSpaceRanger {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getSpaceRangerConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
} ## end sub performSpaceRanger

1;

