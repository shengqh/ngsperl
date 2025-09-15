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
use Storable qw(dclone);
use scRNA::Modules;

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(initializeSpatialTranscriptomeDefaultOptions 
  performSpatialTranscriptome)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeSpatialTranscriptomeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );
  initDefaultValue( $def, "perform_preprocessing", 0);

  return $def;
}

sub getSpatialTranscriptome {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeSpatialTranscriptomeDefaultOptions($def);

  my $project_name = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $tasks = [@$individual, @$summary];

  my $target_dir = getValue($def, "target_dir");

  if ($def->{perform_space_ranger}){
    my $spaceranger_task = "spaceranger";
    my $fastq_folder = getValue($def, "fastq_folder");
    my $spaceranger_transcriptome = getValue($def, "spaceranger_transcriptome");
    my $spaceranger_probe_set = getValue($def, "spaceranger_probe_set");
    my $spaceranger_jobmode = getValue($def, "spaceranger_jobmode", "slurm");
    
    $config->{files} = getValue($def, "files");
    $config->{H_E_brightfield_images} = getValue($def, "H_E_brightfield_images");
    $config->{manual_alignment_json} = $def->{"manual_alignment_json"};

    my $json_option = "";
    if(defined $def->{manual_alignment_json}){
      $json_option = "--loupe-alignment __FILE4__";
    }

    if(defined $def->{id_to_samples}){
      $config->{id_to_samples} = getValue($def, "id_to_samples");
    }else{
      my $all_samples = [];
      for my $k (keys %{$config->{files}}){
        push(@$all_samples, $k);
      }
      $config->{id_to_samples} = { map { $_ => [$_] } @$all_samples };
    }

    $config->{$spaceranger_task} = {
      class => "CQS::ProgramWrapperOneToOne",
      target_dir => "$target_dir/$spaceranger_task",
      perform => 1,
      program => "",
      check_program => 0,
      option => "

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

",
      parameterSampleFile1_ref => "files",
      parameterSampleFile2_ref => "H_E_brightfield_images",
      parameterSampleFile3_ref => "id_to_samples",
      sh_direct => 1,
      no_docker => 1,
      no_output => 1,
      output_ext => "__NAME__/outs/web_summary.html,__NAME__/outs/segmented_outputs/filtered_feature_cell_matrix.h5",
      pbs => {
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "40gb"
      }
    };
    if(defined $def->{manual_alignment_json}){
      $config->{$spaceranger_task}{parameterSampleFile4_ref} = "manual_alignment_json";
    }

    push (@$tasks, $spaceranger_task);

    my $copy_report_task = "${spaceranger_task}_report";
    $config->{$copy_report_task} = {
      class => "CQS::ProgramWrapperOneToOne",
      target_dir => "$target_dir/$copy_report_task",
      perform => 1,
      program => "",
      check_program => 0,
      option => "

cp __FILE__ __NAME__.web_summary.html

",
      parameterSampleFile1_ref => [ $spaceranger_task, "web_summary.html" ],
      sh_direct => 1,
      no_docker => 1,
      no_output => 1,
      output_ext => ".web_summary.html",
      pbs => {
        "nodes"    => "1:ppn=1",
        "walltime" => "1",
        "mem"      => "2gb"
      }
    };
    push (@$tasks, $copy_report_task);
  }

  if ($def->{perform_RCTD}){
    my $rctd_task = "RCTD";

    my $binsize = "8";
    my $RCTD_thread = getValue($def, "RCTD_thread", 16);

    $config->{$rctd_task} = {
      class => "CQS::IndividualR",
      target_dir => "$target_dir/$rctd_task",
      perform => 1,
      option => "",
      rtemplate => "../scRNA/Deconvolution_functions.R,../scRNA/Deconvolution_RCTD.r",
      rReportTemplate => "../scRNA/Deconvolution_RCTD.rmd;reportFunctions.R",
      run_rmd_independent => 1,
      rmd_ext => ".Deconvolution_RCTD.html",
      parameterSampleFile1_ref => "files",
      parameterSampleFile2 => {
        "bin.size" => $binsize,
        "RCTD_thread" => $RCTD_thread,
        "email" => getValue($def, "email"),
        "affiliation" => getValue($def, "affiliation", "CQS/Biostatistics, VUMC"),
      },
      parameterFile1 => getValue($def, "reference"),
      sh_direct => 0,
      no_docker => getValue($def, "no_docker", 0),
      output_ext => ".post_RCTD.RDS,.post_RCTD.valid.RDS",
      pbs => {
        "nodes"    => "1:ppn=${RCTD_thread}",
        "walltime" => "24",
        "mem"      => "40gb"
      }
    };
    push (@$tasks, $rctd_task);

    my $rctd_report_task = $rctd_task . "_report";
    $config->{$rctd_report_task} = {
      class => "CQS::UniqueRmd",
      target_dir => "$target_dir/$rctd_report_task",
      perform => 1,
      option => "",
      report_rmd_file => "../scRNA/Deconvolution_RCTD_report.rmd",
      additional_rmd_files => "../CQS/reportFunctions.R",
      parameterSampleFile1_ref => [ $rctd_task, ".post_RCTD.RDS" ],
      parameterSampleFile2 => {
        email => getValue($def, "email"),
        affiliation => getValue($def, "affiliation", "CQS/Biostatistics, VUMC"),
      },
      output_file_ext => ".RCTD.html",
      output_other_ext => ".RCTD.html",
      sh_direct => 0,
      no_docker => getValue($def, "no_docker", 0),
      pbs => {
        "nodes"    => "1:ppn=1",
        "walltime" => "8",
        "mem"      => "40gb"
      }
    };
    push (@$tasks, $rctd_report_task);  

    my $singlet_task = $rctd_task . "_singlet";
    $config->{$singlet_task} = {
      class => "CQS::IndividualR",
      target_dir => "$target_dir/$singlet_task",
      perform => 1,
      option => "",
      rtemplate => "../scRNA/Deconvolution_RCTD_singlet.r",
      parameterSampleFile1_ref => [ $rctd_task, ".post_RCTD.RDS" ],
      parameterSampleFile2 => {
        "bin.size" => $binsize,
      },
      sh_direct => 0,
      no_docker => getValue($def, "no_docker", 0),
      output_ext => ".post_RCTD.singlet.RDS",
      pbs => {
        "nodes"    => "1:ppn=1",
        "walltime" => "4",
        "mem"      => "20gb"
      }
    };
    push (@$tasks, $singlet_task);
  }

  $config->{sequencetask} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => $tasks,
    },
    sh_direct => 0,
    cluster   => $cluster,
    pbs       => {
      "nodes"     => "1:ppn=" . $def->{max_thread},
      "walltime"  => $def->{sequencetask_run_time},
      "mem"       => "40gb"
    },
  };

  return ($config);
}

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
}

1;

