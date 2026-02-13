#!/usr/bin/perl
package Pipeline::SpatialTranscriptome;

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
use Pipeline::SpaceRanger;

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
  initDefaultValue( $def, "perform_SpaGene",       0 );
  initDefaultValue( $def, "nCount_cutoff",         100 );
  initDefaultValue( $def, "cluster_algorithm",     4 );         #Leiden algorithm

  if ( defined( $def->{bubblemap_width_in} ) ) {
    initDefaultValue( $def, "bubblemap_width", $def->{bubblemap_width_in} * 300 );
  }
  elsif ( defined( $def->{bubblemap_width} ) ) {
    initDefaultValue( $def, "bubblemap_width_in", $def->{bubblemap_width} / 300 );
  }

  initDefaultValue( $def, "MEcell_resolutions", [ 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5 ] );

  return $def;
} ## end sub initializeSpatialTranscriptomeDefaultOptions


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
    my $raw_spaceranger_summary_task     = add_SpaceRanger_summary( $config, $def, $target_dir, "raw_spaceranger", $tasks, $raw_spaceranger_copy_report_task );
  } ## end if ( defined $def->{raw_spaceranger_metrics...})

  if ( $def->{perform_VisiumHD} ) {
    my $VisiumHD_image_features_task = undef;
    if ( $def->{extract_visiumhd_image_features} ) {
      $VisiumHD_image_features_task = "VisiumHD_cell_image_features";

      my $image_features_script = dirname(__FILE__) . "/../scRNA/extract_visiumhd_image_features.py";
      $config->{$VisiumHD_image_features_task} = {
        class         => "CQS::ProgramWrapperOneToOne",
        perform       => 1,
        target_dir    => "${target_dir}/$VisiumHD_image_features_task",
        program       => "",
        check_program => 0,
        option        => "
python3 $image_features_script \\
  --cell_geojson '__FILE__' \\
  --nucleus_geojson '__FILE2__' \\
  --scales '__FILE3__' \\
  --image '__FILE4__' \\
  --output '__NAME__.resnet_features.parquet' 
",
        parameterSampleFile1 => getValue($def, "cell_geojson_files"),
        parameterSampleFile2 => getValue($def, "nucleus_geojson_files"),
        parameterSampleFile3 => getValue($def, "scales_files"),
        parameterSampleFile4 => getValue($def, "image_files"),
        output_ext               => ".resnet_features.parquet",
        docker_prefix            => "visiumhd_",
        no_output                => 1,
        sh_direct                => 0,
        pbs                      => {
          "nodes"    => "1:ppn=4",
          "walltime" => "10:00:00",
          "mem"      => "80gb"
        },
      };
      push( @$tasks, $VisiumHD_image_features_task );
    } ## end if ( $def->{extract_visiumhd_image_features...})

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
        "email"           => getValue( $def, "email" ),
        "affiliation"     => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
        "Mtpattern"       => "^MT-|^Mt-|^mt-",
        "bin_size"        => $bin_size,
        "nCount_cutoff"   => getValue( $def, "nCount_cutoff" ),
        "nFeature_cutoff" => getValue( $def, "nFeature_cutoff" ),
        "mt_cutoff"       => getValue( $def, "mt_cutoff" ),
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

    my $filter_task = "VisiumHD_filter_cellonly";
    $config->{$filter_task} = {
      class           => "CQS::IndividualRmd",
      target_dir      => "$target_dir/$filter_task",
      perform         => 1,
      rReportTemplate => "../scRNA/VisiumHD_filter.Rmd;reportFunctions.R;../scRNA/scRNA_func.r;../scRNA/spatial_plotting_functions.R",
      option          => "

Rscript --vanilla  -e \"library('rmarkdown');rmarkdown::render('VisiumHD_filter.Rmd', output_file='__OUTPUT__')\"

",
      parameterSampleFile1_ref => [ $qc_task, ".rds" ],
      parameterSampleFile2     => {
        "email"           => getValue( $def, "email" ),
        "affiliation"     => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
        "nCount_cutoff"   => getValue( $def, "nCount_cutoff" ),
        "nFeature_cutoff" => getValue( $def, "nFeature_cutoff" ),
        "mt_cutoff"       => getValue( $def, "mt_cutoff" ),
        "bin_size"        => "Spatial.Polygons",
      },
      no_prefix             => 1,
      sh_direct             => 0,
      no_docker             => getValue( $def, "no_docker", 0 ),
      output_to_same_folder => 0,
      output_file_prefix    => ".filtered.html",
      output_ext            => ".filtered.html",
      output_other_ext      => ".filtered.rds",
      pbs                   => {
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "40gb"
      }
    };
    push( @$tasks, $filter_task );

    my $spagene_task = undef;
    if ( $def->{perform_SpaGene} ) {
      $spagene_task = "SpaGene";
      $config->{$spagene_task} = {
        class                    => "CQS::IndividualR",
        target_dir               => "$target_dir/$spagene_task",
        perform                  => 1,
        option                   => "",
        rtemplate                => "reportFunctions.R,../scRNA/SpaGene.r",
        parameterSampleFile1_ref => [ $filter_task, ".rds" ],
        parameterSampleFile2     => {
          "email"         => getValue( $def, "email" ),
          "affiliation"   => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
          "nCount_cutoff" => getValue( $def, "nCount_cutoff" ),
          "LRpair_file"   => getValue( $def, "SpaGene_LRpair_file" ),
        },
        sh_direct  => 0,
        no_docker  => getValue( $def, "SpaGene_no_docker", 1 ),
        output_ext => ".spa.rds,.spa_pattern.rds,.spa_lr.rds",
        pbs        => {
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "300gb"
        }
      };
      push( @$tasks, $spagene_task );
    } ## end if ( $def->{perform_SpaGene...})

    my $rctd_polygons_task = undef;
    if ( $def->{perform_RCTD} ) {
      my $rctd_tasks = {};
      # We cannot set RCTD too many thread. It is very easy to get error in cluster.
      my $RCTD_thread = getValue( $def, "RCTD_thread", 8 );
      #my $assays      = [ 'Spatial.Polygons', 'Spatial.008um' ];
      # we will focus on Polygons only
      my $assays = ['Spatial.Polygons'];
      for my $assay ( @{$assays} ) {
        my $rctd_task = "RCTD_$assay";
        if ( $assay eq 'Spatial.Polygons' ) {
          $rctd_polygons_task = $rctd_task;
        }
        $config->{$rctd_task} = {
          class                    => "CQS::IndividualR",
          target_dir               => "$target_dir/$rctd_task",
          perform                  => 1,
          option                   => "",
          rtemplate                => "reportFunctions.R,../scRNA/Deconvolution_functions.R,../scRNA/Deconvolution_RCTD_obj.r",
          parameterSampleFile1_ref => [ $filter_task, ".rds" ],
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

      my $vis_task = "RCTD_visualization";
      $config->{$vis_task} = {
        class                    => "CQS::IndividualRmd",
        target_dir               => "$target_dir/$vis_task",
        perform                  => 1,
        option                   => "",
        rReportTemplate          => "../scRNA/Deconvolution_RCTD_visualization.rmd,../scRNA/scRNA_func.r;../CQS/reportFunctions.R;../scRNA/Deconvolution_functions.R",
        parameterSampleFile1_ref => [ $rctd_tasks->{'Spatial.Polygons'}, ".RCTD.obj.rds" ],
        parameterSampleFile2     => {
          email       => getValue( $def, "email" ),
          affiliation => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
        },
        output_ext => ".RCTD_visualization.html",
        sh_direct  => 0,
        no_docker  => getValue( $def, "no_docker", 0 ),
        pbs        => {
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "40gb"
        }
      };
      if ( $rctd_tasks->{'Spatial.008um'} ) {
        $config->{$vis_task}->{parameterSampleFile3_ref} = [ $rctd_tasks->{'Spatial.008um'}, ".RCTD.obj.rds" ];
      }
      push( @$tasks, $vis_task );

      my $vis_summary_task = "${vis_task}_summary";
      $config->{$vis_summary_task} = {
        class                    => "CQS::UniqueRmd",
        target_dir               => "$target_dir/$vis_summary_task",
        perform                  => 1,
        option                   => "",
        report_rmd_file          => "../scRNA/Deconvolution_RCTD_visualization_summary.Rmd",
        additional_rmd_files     => "../CQS/reportFunctions.R",
        parameterSampleFile1_ref => [ $vis_task, ".html" ],
        parameterSampleFile2     => {
          task_name   => getValue( $def, "task_name" ),
          email       => getValue( $def, "email" ),
          affiliation => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
        },
        output_file_ext  => ".RCTD_visualization_summary.html",
        output_other_ext => ".RCTD_visualization_summary.html",
        sh_direct        => 0,
        no_docker        => getValue( $def, "no_docker", 0 ),
        pbs              => {
          "nodes"    => "1:ppn=1",
          "walltime" => "8",
          "mem"      => "40gb"
        }
      };
      push( @$tasks, $vis_summary_task );

    } ## end if ( $def->{perform_RCTD...})

    my $azimuth_task = undef;
    if ( getValue( $def, "perform_Azimuth", 0 ) ) {
      $azimuth_task = "Azimuth";
      $config->{$azimuth_task} = {
        class                    => "CQS::IndividualR",
        target_dir               => "$target_dir/$azimuth_task",
        perform                  => 1,
        option                   => "",
        rtemplate                => "../CQS/reportFunctions.R,../scRNA/scRNA_func.r,../scRNA/spatial_azimuth.r",
        parameterSampleFile1_ref => [ $filter_task, ".rds" ],
        parameterSampleFile2     => {
          task_name       => getValue( $def, "task_name" ),
          email           => getValue( $def, "email" ),
          affiliation     => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
          assay           => "Spatial.Polygons",
          Azimuth_ref     => getValue( $def, "Azimuth_ref" ),
          bubblemap_file  => getValue( $def, "bubblemap_file" ),
          bubblemap_width => getValue( $def, "bubblemap_width" ),
          nCount_cutoff   => getValue( $def, "nCount_cutoff" ),
        },
        sh_direct     => 0,
        docker_prefix => "azimuth_",
        no_docker     => 0,
        output_ext    => ".meta.rds",
        pbs           => {
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "40gb"
        }
      };
      push( @$tasks, $azimuth_task );
    } ## end if ( getValue( $def, "perform_Azimuth"...))

    my $singleR_task = undef;
    if ( getValue( $def, "perform_SingleR", 0 ) ) {
      $singleR_task = "SingleR";
      $config->{$singleR_task} = {
        class                    => "CQS::IndividualR",
        target_dir               => "$target_dir/$singleR_task",
        perform                  => 1,
        option                   => "",
        rtemplate                => "../CQS/reportFunctions.R,../scRNA/scRNA_func.r,../scRNA/SingleR.r",
        parameterSampleFile1_ref => [ $filter_task, ".rds" ],
        parameterSampleFile2     => {
          task_name       => getValue( $def, "task_name" ),
          email           => getValue( $def, "email" ),
          affiliation     => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
          assay           => "Spatial.Polygons",
          bubblemap_file  => getValue( $def, "bubblemap_file" ),
          bubblemap_width => getValue( $def, "bubblemap_width" ),
          nCount_cutoff   => getValue( $def, "nCount_cutoff" ),
          species         => getValue( $def, "species" ),
        },
        sh_direct     => 0,
        docker_prefix => "singleR_",
        no_docker     => 0,
        output_ext    => ".meta.rds",
        pbs           => {
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "40gb"
        }
      };
      push( @$tasks, $singleR_task );
    } ## end if ( getValue( $def, "perform_SingleR"...))

    my $signacx_task = undef;
    if ( getValue( $def, "perform_SignacX", 0 ) ) {
      $signacx_task = "SignacX";
      $config->{$signacx_task} = {
        class                    => "CQS::IndividualR",
        target_dir               => "$target_dir/$signacx_task",
        perform                  => 1,
        option                   => "",
        rtemplate                => "../CQS/reportFunctions.R,../scRNA/scRNA_func.r,../scRNA/SignacX_only.r",
        parameterSampleFile1_ref => [ $filter_task, ".rds" ],
        parameterSampleFile2     => {
          task_name       => getValue( $def, "task_name" ),
          email           => getValue( $def, "email" ),
          affiliation     => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
          assay           => "Spatial.Polygons",
          bubblemap_file  => getValue( $def, "bubblemap_file" ),
          bubblemap_width => getValue( $def, "bubblemap_width" ),
          nCount_cutoff   => getValue( $def, "nCount_cutoff" ),
          species         => getValue( $def, "species" ),
          pca_dims        => getValue( $def, "pca_dims", 30 ),
          reduction       => "pca",
          by_sctransform  => 0,
        },
        sh_direct     => 0,
        docker_prefix => "signacX_",
        no_docker     => 0,
        output_ext    => ".SignacX.rds",
        pbs           => {
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "40gb"
        }
      };
      push( @$tasks, $signacx_task );
    } ## end if ( getValue( $def, "perform_SignacX"...))

    if ( getValue( $def, "perform_MEcell", 0 ) ) {
      my $MEcell_task = "MEcell";
      $config->{$MEcell_task} = {
        class                    => "CQS::IndividualR",
        target_dir               => "$target_dir/$MEcell_task",
        perform                  => 1,
        option                   => "",
        rtemplate                => "../scRNA/scRNA_func.r,../scRNA/spatial_MEcell.r",
        parameterSampleFile1_ref => [ $filter_task, ".rds" ],
        parameterSampleFile2     => {
          task_name     => getValue( $def, "task_name" ),
          email         => getValue( $def, "email" ),
          affiliation   => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
          assay         => "Spatial.Polygons",
          nCount_cutoff => getValue( $def, "nCount_cutoff" ),
        },
        sh_direct        => 0,
        no_docker        => 0,
        output_file_ext  => ".MEcell.rds",
        output_other_ext => ".MEcell.mtx.gz",
        pbs              => {
          "nodes"    => "1:ppn=1",
          "walltime" => "48",
          "mem"      => "40gb"
        }
      };
      push( @$tasks, $MEcell_task );

      my $MEcell_umap_task = "MEcell_umap";
      $config->{$MEcell_umap_task} = {
        class                    => "CQS::IndividualR",
        target_dir               => "$target_dir/$MEcell_umap_task",
        perform                  => 1,
        option                   => "",
        rtemplate                => "../scRNA/scRNA_func.r,../scRNA/spatial_MEcell_umap_graph.r",
        parameterSampleFile1_ref => [ $MEcell_task, ".rds" ],
        parameterSampleFile2     => {
          task_name   => getValue( $def, "task_name" ),
          email       => getValue( $def, "email" ),
          affiliation => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
          assay       => "Spatial.Polygons",
          min_neighbors => getValue( $def, "MEcell_umap_min_neighbors" ),
        },
        sh_direct       => 0,
        no_docker       => 0,
        output_file_ext => ".MEcell_umap.rds",
        pbs             => {
          "nodes"    => "1:ppn=1",
          "walltime" => "48",
          "mem"      => "40gb"
        }
      };
      push( @$tasks, $MEcell_umap_task );

      my $MEcell_resolutions = getValue( $def, "MEcell_resolutions" );
      if ( !is_array($MEcell_resolutions) ) {
        stop("MEcell_resolutions should be array of float.");
      }
      $MEcell_resolutions = [ map { $_ + 0 } @{$MEcell_resolutions} ];
      my $MEcell_resolutions_label = join( "_", @{$MEcell_resolutions} );
      my $assay                    = "Spatial.Polygons";
      my $cluster_algorithm        = getValue( $def, "cluster_algorithm", 4 );
      my $cluster_algorithm_name;
      if ( $cluster_algorithm == 1 ) {
        $cluster_algorithm_name = 'Louvain';
      }
      elsif ( $cluster_algorithm == 2 ) {
        $cluster_algorithm_name = 'LouvainRefine';
      }
      elsif ( $cluster_algorithm == 3 ) {
        $cluster_algorithm_name = 'SLM';
      }
      elsif ( $cluster_algorithm == 4 ) {
        $cluster_algorithm_name = 'Leiden';
      }
      else {
        die "Unsupported cluster_algorithm $cluster_algorithm. Supported values are 1 (Louvain), 2 (LouvainRefine), 3 (SLM), 4 (Leiden).";
      }

      my $MEcell_cluster_tasks = [];
      for ( my $i = 0; $i < scalar( @{$MEcell_resolutions} ); $i++ ) {
        my $cur_resolution = $MEcell_resolutions->[$i];

        my $MEcell_cluster_task = "${MEcell_task}_cluster_${cluster_algorithm_name}_res.${cur_resolution}";
        $config->{$MEcell_cluster_task} = {
          class                    => "CQS::IndividualR",
          target_dir               => "$target_dir/$MEcell_cluster_task",
          perform                  => 1,
          rtemplate                => "reportFunctions.R;../scRNA/scRNA_func.r;../scRNA/spatial_MEcell_one_cluster.r",
          option                   => "",
          parameterSampleFile1_ref => [ $MEcell_task, ".rds" ],
          parameterSampleFile2     => {
            email                  => getValue( $def, "email" ),
            affiliation            => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
            assay                  => $assay,
            markers_file           => getValue( $def, "markers_file" ),
            curated_markers_file   => getValue( $def, "curated_markers_file" ),
            summary_layer_file     => getValue( $def, "summary_layer_file" ),
            remove_subtype         => getValue( $def, "remove_subtype" ),
            HLA_panglao5_file      => getValue( $def, "HLA_panglao5_file" ),
            bubblemap_file         => getValue( $def, "bubblemap_file" ),
            bubblemap_width_in     => getValue( $def, "bubblemap_width_in", 8 ),
            MEcell_resolution      => $cur_resolution,
            species                => getValue( $def, "species" ),
            cluster_algorithm      => $cluster_algorithm,
            cluster_algorithm_name => $cluster_algorithm_name,
          },
          parameterSampleFile3_ref => $azimuth_task,
          parameterSampleFile4_ref => $rctd_polygons_task,
          no_prefix                => 1,
          sh_direct                => 0,
          no_docker                => getValue( $def, "no_docker", 0 ),
          output_to_same_folder    => 0,
          output_ext               => ".$assay.MEcell.res.$cur_resolution.$cluster_algorithm_name.meta.rds",
          #output_other_ext         => ".MEcell_cluster.html",
          pbs => {
            "nodes"    => "1:ppn=1",
            "walltime" => "48",
            "mem"      => "40gb"
          }
        };
        push( @$tasks,                $MEcell_cluster_task );
        push( @$MEcell_cluster_tasks, $MEcell_cluster_task );
      } ## end for ( my $i = 0; $i < scalar...)

      my $MEcell_cluster_report_task = "${MEcell_task}_cluster_${cluster_algorithm_name}_report";
      $config->{$MEcell_cluster_report_task} = {
        class                    => "CQS::IndividualRmd",
        target_dir               => "$target_dir/$MEcell_cluster_report_task",
        perform                  => 1,
        rReportTemplate          => "../scRNA/spatial_MEcell_cluster.rmd;reportFunctions.R;../scRNA/scRNA_func.r",
        option                   => "",
        parameterSampleFile1_ref => [ $MEcell_task, ".rds" ],
        parameterSampleFile2     => {
          email                  => getValue( $def, "email" ),
          affiliation            => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
          assay                  => $assay,
          markers_file           => getValue( $def, "markers_file" ),
          curated_markers_file   => getValue( $def, "curated_markers_file" ),
          summary_layer_file     => getValue( $def, "summary_layer_file" ),
          remove_subtype         => getValue( $def, "remove_subtype" ),
          HLA_panglao5_file      => getValue( $def, "HLA_panglao5_file" ),
          bubblemap_file         => getValue( $def, "bubblemap_file" ),
          bubblemap_width_in     => getValue( $def, "bubblemap_width_in", 8 ),
          MEcell_resolutions     => $MEcell_resolutions,
          species                => getValue( $def, "species" ),
          cluster_algorithm      => $cluster_algorithm,
          cluster_algorithm_name => $cluster_algorithm_name,
        },
        parameterSampleFile3_ref => $MEcell_cluster_tasks,
        parameterSampleFile4_ref => $azimuth_task,
        parameterSampleFile5_ref => $rctd_polygons_task,
        parameterSampleFile6_ref => $singleR_task,
        parameterSampleFile7_ref => $signacx_task,
        parameterSampleFile8_ref => $MEcell_umap_task,
        parameterSampleFile9_ref => $VisiumHD_image_features_task,
        no_prefix                => 1,
        sh_direct                => 0,
        no_docker                => getValue( $def, "no_docker", 0 ),
        output_to_same_folder    => 0,
        output_ext               => ".MEcell_cluster.html",
        pbs                      => {
          "nodes"    => "1:ppn=1",
          "walltime" => "48",
          "mem"      => "40gb"
        }
      };
      push( @$tasks, $MEcell_cluster_report_task );

      my $MEcell_cluster_summary_task = "${MEcell_cluster_report_task}_summary";
      $config->{$MEcell_cluster_summary_task} = {
        class                    => "CQS::UniqueRmd",
        target_dir               => "$target_dir/$MEcell_cluster_summary_task",
        perform                  => 1,
        option                   => "",
        report_rmd_file          => "../scRNA/spatial_MEcell_cluster_summary.Rmd",
        additional_rmd_files     => "../CQS/reportFunctions.R",
        parameterSampleFile1_ref => [ $MEcell_cluster_report_task, ".html" ],
        parameterSampleFile2     => {
          task_name   => getValue( $def, "task_name" ),
          email       => getValue( $def, "email" ),
          affiliation => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
        },
        output_file_ext  => ".MEcell_cluster_summary.html",
        output_other_ext => ".MEcell_cluster_summary.html",
        sh_direct        => 0,
        no_docker        => getValue( $def, "no_docker", 0 ),
        pbs              => {
          "nodes"    => "1:ppn=1",
          "walltime" => "8",
          "mem"      => "40gb"
        }
      };
      push( @$tasks, $MEcell_cluster_summary_task );

    } ## end if ( getValue( $def, "perform_MEcell"...))

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

