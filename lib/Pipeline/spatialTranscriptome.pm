#!/usr/bin/perl
package Pipeline::spatialTranscriptome;

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
    qw(initializeSpatialTranscriptomeDefaultOptions
        performSpatialTranscriptome)
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';


sub initializeSpatialTranscriptomeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType",             "FAIL" );
  initDefaultValue( $def, "cluster",               "slurm" );
  initDefaultValue( $def, "perform_preprocessing", 0 );
  initDefaultValue( $def, "perform_VisiumHD",      1 );
  initDefaultValue( $def, "perform_RCTD",          1 );

  return $def;
} ## end sub initializeSpatialTranscriptomeDefaultOptions


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


sub add_spaceranger_summary {
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
} ## end sub add_spaceranger_summary


sub getSpatialTranscriptome {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeSpatialTranscriptomeDefaultOptions($def);

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
    my $raw_spaceranger_summary_task     = add_spaceranger_summary( $config, $def, $target_dir, "raw_spaceranger", $tasks, $raw_spaceranger_copy_report_task );
  } ## end if ( defined $def->{raw_spaceranger_metrics...})

  if ( $def->{perform_space_ranger} ) {
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
    my $spaceranger_summary_task     = add_spaceranger_summary( $config, $def, $target_dir, $spaceranger_task, $tasks, $spaceranger_copy_report_task );

    # spaceranger should not be performed together with VisiumHD downstream analysis
    $def->{perform_VisiumHD} = 0;
    $def->{perform_RCTD}     = 0;
  } ## end if ( $def->{perform_space_ranger...})

  if ( $def->{perform_VisiumHD} ) {
    my $bin_size       = getValue( $def, "bin_size", "8,polygons" );
    my $bin_size_label = $bin_size;
    $bin_size_label =~ s/,/_/;

    my $qc_task = "VisiumHD_qc";
    $config->{$qc_task} = {
      class           => "CQS::IndividualRmd",
      target_dir      => "$target_dir/$qc_task",
      perform         => 1,
      rReportTemplate => "../scRNA/VisiumHD_qc.Rmd;reportFunctions.R;../scRNA/scRNA_func.r;../scRNA/spatial_plotting_functions.R",
      option          => "

Rscript --vanilla  -e \"library('rmarkdown');rmarkdown::render('VisiumHD_qc.Rmd', output_file='__OUTPUT__')\"

",
      parameterSampleFile1_ref => "files",
      parameterSampleFile2     => {
        "email"       => getValue( $def, "email" ),
        "affiliation" => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
        "Mtpattern"   => "^MT-|^Mt-|^mt-",
        "bin_size"    => $bin_size,
      },
      no_prefix             => 1,
      sh_direct             => 0,
      no_docker             => getValue( $def, "no_docker", 0 ),
      output_to_same_folder => 0,
      output_file_prefix    => ".$bin_size_label.html",
      output_ext            => ".$bin_size_label.html",
      output_other_ext      => ".$bin_size_label.rds",
      pbs                   => {
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "40gb"
      }
    };
    push( @$tasks, $qc_task );

    my $qc_summary_task = "${qc_task}_summary";
    $config->{$qc_summary_task} = {
      class                    => "CQS::UniqueRmd",
      target_dir               => "$target_dir/$qc_summary_task",
      perform                  => 1,
      option                   => "",
      report_rmd_file          => "../scRNA/VisiumHD_qc_summary.Rmd",
      additional_rmd_files     => "../CQS/reportFunctions.R",
      parameterSampleFile1_ref => [ $qc_task, ".html" ],
      parameterSampleFile2     => {
        task_name   => getValue( $def, "task_name" ),
        email       => getValue( $def, "email" ),
        affiliation => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
      },
      output_file_ext  => ".VisiumHD_qc_summary.html",
      output_other_ext => ".VisiumHD_qc_summary.html",
      sh_direct        => 0,
      no_docker        => getValue( $def, "no_docker", 0 ),
      pbs              => {
        "nodes"    => "1:ppn=1",
        "walltime" => "8",
        "mem"      => "40gb"
      }
    };
    push( @$tasks, $qc_summary_task );

    if ( $def->{perform_RCTD} ) {
      # We cannot set RCTD too many thread. It is very easy to get error in cluster.
      my $RCTD_thread = getValue( $def, "RCTD_thread", 8 );
      my $assays      = [ 'Spatial.Polygons', 'Spatial.008um' ];
      my $rctd_tasks  = {};
      for my $assay ( @{$assays} ) {
        my $rctd_task = "RCTD_$assay";
        $config->{$rctd_task} = {
          class                    => "CQS::IndividualR",
          target_dir               => "$target_dir/$rctd_task",
          perform                  => 1,
          option                   => "",
          rtemplate                => "reportFunctions.R,../scRNA/Deconvolution_functions.R,../scRNA/Deconvolution_RCTD_obj.r",
          parameterSampleFile1_ref => [ $qc_task, ".rds" ],
          parameterSampleFile2     => {
            "assay"       => $assay,
            "RCTD_thread" => $RCTD_thread,
            "email"       => getValue( $def, "email" ),
            "affiliation" => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
          },
          parameterFile1 => getValue( $def, "reference" ),
          sh_direct      => 0,
          no_docker      => getValue( $def, "no_docker", 0 ),
          output_ext     => ".$assay.RCTD.obj.rds",
          pbs            => {
            "nodes"    => "1:ppn=${RCTD_thread}",
            "walltime" => "24",
            "mem"      => "80gb"
          }
        };
        push( @$tasks, $rctd_task );
        $rctd_tasks->{$assay} = $rctd_task;
      } ## end for my $assay ( @{$assays...})

      my $comparison_task = "RCTD_comparison";
      $config->{$comparison_task} = {
        class                    => "CQS::IndividualRmd",
        target_dir               => "$target_dir/$comparison_task",
        perform                  => 1,
        option                   => "",
        rReportTemplate          => "../scRNA/Deconvolution_RCTD_comparison.rmd,../scRNA/scRNA_func.r;../CQS/reportFunctions.R;../scRNA/Deconvolution_functions.R",
        parameterSampleFile1_ref => [ $rctd_tasks->{'Spatial.Polygons'}, ".RCTD.obj.rds" ],
        parameterSampleFile2     => {
          email       => getValue( $def, "email" ),
          affiliation => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
        },
        parameterSampleFile3_ref => [ $rctd_tasks->{'Spatial.008um'}, ".RCTD.obj.rds" ],
        output_file_ext          => ".RCTD.comparison.html",
        output_other_ext         => ".RCTD.comparison.html",
        sh_direct                => 0,
        no_docker                => getValue( $def, "no_docker", 0 ),
        pbs                      => {
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "40gb"
        }
      };
      push( @$tasks, $comparison_task );

    } ## end if ( $def->{perform_RCTD...})
  } ## end if ( $def->{perform_VisiumHD...})

  #   if ( $def->{perform_segment_report} ) {
  #     my $segment_report_task = "segment_report";
  #     $config->{$segment_report_task} = {
  #       class           => "CQS::IndividualRmd",
  #       target_dir      => "$target_dir/$segment_report_task",
  #       perform         => 1,
  #       program         => "",
  #       check_program   => 0,
  #       rReportTemplate => "../scRNA/spaceranger_cell_individual.Rmd;reportFunctions.R;../scRNA/scRNA_func.r;../scRNA/spaceranger_cell_individual_func.r",
  #       option          => "

  # Rscript --vanilla  -e \"library('rmarkdown');rmarkdown::render('spaceranger_cell_individual.Rmd', output_file='__OUTPUT__')\"

  # ",
  #       parameterSampleFile1_ref => $source_ref,
  #       parameterSampleFile2     => {
  #         email                  => getValue( $def, "email" ),
  #         affiliation            => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
  #         task_name              => getValue( $def, "task_name" ),
  #         Mtpattern              => getValue( $def, "Mtpattern" ),
  #         rRNApattern            => getValue( $def, "rRNApattern" ),
  #         regress_by_percent_mt  => getValue( $def, "regress_by_percent_mt" ),
  #         nFeature_cutoff_min    => getValue( $def, "nFeature_cutoff_min" ),
  #         nFeature_cutoff_max    => getValue( $def, "nFeature_cutoff_max" ),
  #         nCount_cutoff          => getValue( $def, "nCount_cutoff" ),
  #         mt_cutoff              => getValue( $def, "mt_cutoff" ),
  #         species                => getValue( $def, "species" ),
  #         resolution             => getValue( $def, "resolution" ),
  #         pca_dims               => getValue( $def, "pca_dims" ),
  #         by_sctransform         => getValue( $def, "by_sctransform" ),
  #         use_sctransform_v2     => getValue( $def, "use_sctransform_v2", 1 ),
  #         markers_file           => getValue( $def, "markers_file" ),
  #         curated_markers_file   => getValue( $def, "curated_markers_file" ),
  #         remove_subtype         => getValue( $def, "remove_subtype" ),
  #         HLA_panglao5_file      => getValue( $def, "HLA_panglao5_file" ),
  #         bubblemap_file         => getValue( $def, "bubblemap_file" ),
  #         celltype_predictmethod => getValue( $def, "celltype_predictmethod", "cta" ),
  #       },
  #       no_prefix             => 1,
  #       sh_direct             => 1,
  #       no_docker             => 1,
  #       output_to_same_folder => 0,
  #       output_file_prefix    => ".segment_report.html",
  #       output_ext            => ".segment_report.html",
  #       pbs                   => {
  #         "nodes"    => "1:ppn=1",
  #         "walltime" => "24",
  #         "mem"      => "40gb"
  #       }
  #     };
  #     push( @$tasks, $segment_report_task );
  #   } ## end if ( $def->{perform_segment_report...})

  # if ( $def->{perform_RCTD} ) {
  #   my $rctd_task = "RCTD";

  #   my $binsize     = "8";
  #   my $RCTD_thread = getValue( $def, "RCTD_thread", 16 );

  #   $config->{$rctd_task} = {
  #     class                    => "CQS::IndividualR",
  #     target_dir               => "$target_dir/$rctd_task",
  #     perform                  => 1,
  #     option                   => "",
  #     rtemplate                => "../scRNA/Deconvolution_functions.R,../scRNA/Deconvolution_RCTD.r",
  #     rReportTemplate          => "../scRNA/Deconvolution_RCTD.rmd;reportFunctions.R",
  #     run_rmd_independent      => 1,
  #     rmd_ext                  => ".Deconvolution_RCTD.html",
  #     parameterSampleFile1_ref => "files",
  #     parameterSampleFile2     => {
  #       "bin.size"    => $binsize,
  #       "RCTD_thread" => $RCTD_thread,
  #       "email"       => getValue( $def, "email" ),
  #       "affiliation" => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
  #     },
  #     parameterFile1 => getValue( $def, "reference" ),
  #     sh_direct      => 0,
  #     no_docker      => getValue( $def, "no_docker", 0 ),
  #     output_ext     => ".post_RCTD.RDS,.post_RCTD.valid.RDS",
  #     pbs            => {
  #       "nodes"    => "1:ppn=${RCTD_thread}",
  #       "walltime" => "24",
  #       "mem"      => "40gb"
  #     }
  #   };
  #   push( @$tasks, $rctd_task );

  #   my $rctd_report_task = $rctd_task . "_report";
  #   $config->{$rctd_report_task} = {
  #     class                    => "CQS::UniqueRmd",
  #     target_dir               => "$target_dir/$rctd_report_task",
  #     perform                  => 1,
  #     option                   => "",
  #     report_rmd_file          => "../scRNA/Deconvolution_RCTD_report.rmd",
  #     additional_rmd_files     => "../CQS/reportFunctions.R",
  #     parameterSampleFile1_ref => [ $rctd_task, ".post_RCTD.RDS" ],
  #     parameterSampleFile2     => {
  #       task_name   => getValue( $def, "task_name" ),
  #       email       => getValue( $def, "email" ),
  #       affiliation => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
  #     },
  #     output_file_ext  => ".RCTD.html",
  #     output_other_ext => ".RCTD.html",
  #     sh_direct        => 0,
  #     no_docker        => getValue( $def, "no_docker", 0 ),
  #     pbs              => {
  #       "nodes"    => "1:ppn=1",
  #       "walltime" => "8",
  #       "mem"      => "40gb"
  #     }
  #   };
  #   push( @$tasks, $rctd_report_task );

  #   my $singlet_task = $rctd_task . "_singlet";
  #   $config->{$singlet_task} = {
  #     class                    => "CQS::IndividualR",
  #     target_dir               => "$target_dir/$singlet_task",
  #     perform                  => 1,
  #     option                   => "",
  #     rtemplate                => "../scRNA/Deconvolution_RCTD_singlet.r",
  #     parameterSampleFile1_ref => [ $rctd_task, ".post_RCTD.RDS" ],
  #     parameterSampleFile2     => { "bin.size" => $binsize, },
  #     sh_direct                => 0,
  #     no_docker                => getValue( $def, "no_docker", 0 ),
  #     output_ext               => ".post_RCTD.singlet.RDS",
  #     pbs                      => {
  #       "nodes"    => "1:ppn=1",
  #       "walltime" => "4",
  #       "mem"      => "20gb"
  #     }
  #   };
  #   push( @$tasks, $singlet_task );
  # } ## end if ( $def->{perform_RCTD...})

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
} ## end sub getSpatialTranscriptome


sub performSpatialTranscriptome {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getSpatialTranscriptome($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
} ## end sub performSpatialTranscriptome

1;

