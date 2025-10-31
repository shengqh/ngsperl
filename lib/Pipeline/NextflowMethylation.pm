#!/usr/bin/perl
package Pipeline::NextflowMethylation;

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
use Pipeline::MethylationUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performNextflowMethylation)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';


sub initializeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType",             "FAIL" );
  initDefaultValue( $def, "cluster",               "slurm" );
  initDefaultValue( $def, "perform_preprocessing", 0 );
  initDefaultValue( $def, "aligner",               "bismark" );

  initDefaultValue( $def, "methylation_mincov", "20" );

  return $def;
} ## end sub initializeDefaultOptions


sub getNextflowMethylationConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

  my $task_name = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $tasks = [ @$individual, @$summary ];

  my $target_dir = $def->{target_dir};

  my $convered_bed     = getValue( $def, "convered_bed" );
  my $igenomes_base    = getValue( $def, "igenomes_base" );
  my $genome           = getValue( $def, "genome" );
  my $aligner          = getValue( $def, "aligner" );
  my $nextflow_config  = getValue( $def, "nextflow_config" );
  my $nextflow_main_nf = getValue( $def, "nextflow_main_nf" );
  my $sh_direct        = getValue( $def, "sh_direct" );

  my $nextflow_run_mode = getValue( $def, "nextflow_run_mode", "local" );

  my $nextflow_methylseq_task = "nextflow_methylseq";

  my $perform_downstream_analysis = getValue( $def, "perform_downstream_analysis", 0 );

  my $methylkitprep_task = "MethylKitPreparation";

  if ( $nextflow_run_mode eq "slurm" ) {
    # in slurm mode, we can let nextflow to handle all samples together, the we can only need to run one main script at command line.
    $config->{$nextflow_methylseq_task} = {
      class      => "CQS::ProgramWrapper",
      target_dir => $target_dir . "/$nextflow_methylseq_task",
      option     => "
nextflow run $nextflow_main_nf \\
  -config $nextflow_config \\
  -profile singularity \\
  --input fileList1.list.csv \\
  --outdir . \\
  --genome $genome \\
  --igenomes_base $igenomes_base \\
  --aligner $aligner \\
  --run_targeted_sequencing \\
  --target_regions_file $convered_bed \\
  --collecthsmetrics

status=\$?
if [ \$status -ne 0 ]; then
  echo \"Error: nextflow run nf-core/methylseq failed with status \$status\"
  exit \$status
else
  echo \"nextflow run nf-core/methylseq completed successfully\"
  rm -rf work
fi

# for list input files: ",
      program                             => "",
      check_program                       => 0,
      parameterSampleFile1_ref            => $source_ref,
      parameterSampleFile1_header         => "sample,fastq_1,fastq_2,genome",
      parameterSampleFile1_join_delimiter => ",",
      parameterSampleFile1_fileFirst      => 0,
      parameterSampleFile1_col_delimiter  => ",",
      parameterSampleFile1_fileSuffix     => ".csv",
      parameterSampleFile1_suffix         => ",",
      no_prefix                           => 1,
      no_output                           => 1,
      output_ext                          => "bismark/summary/bismark_summary_report.html",
      samplename_in_result                => 0,
      output_to_same_folder               => 1,
      sh_direct                           => $sh_direct,
      no_docker                           => 1,
      pbs                                 => {
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "10gb",
      },
    };

    push @$tasks, $nextflow_methylseq_task;

    if ($perform_downstream_analysis) {
      my $nextflow_methylseq_bismark_dir = $config->{"nextflow_methylseq"}{target_dir} . "/result/bismark/methylation_calls/methylation_coverage";
      $config->{$methylkitprep_task} = {
        class      => "CQS::ProgramWrapperOneToOne",
        target_dir => $target_dir . "/$methylkitprep_task",
        option     => "
rm -f __NAME__.bismark.cov.gz

ln -s $nextflow_methylseq_bismark_dir/__NAME___1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz __NAME__.bismark.cov.gz

# __FILE__
",
        program                  => "",
        check_program            => 0,
        parameterSampleFile1_ref => $source_ref,
        parameterSampleFile2_ref => $nextflow_methylseq_task, # just to make dependency
        no_prefix                => 1,
        no_output                => 1,
        output_ext               => ".bismark.cov.gz",
        output_to_same_folder    => 1,
        sh_direct                => 1,
        pbs                      => {
          "nodes"    => "1:ppn=1",
          "walltime" => "1",
          "mem"      => "2gb",
        },
      };

      push( @$tasks, $methylkitprep_task );
    } ## end if ($perform_downstream_analysis)
  } ## end if ( $nextflow_run_mode...)
  else {
    # in local mode, each sample runs indpendently
    $config->{$nextflow_methylseq_task} = {
      class      => "CQS::ProgramWrapperOneToOne",
      target_dir => $target_dir . "/$nextflow_methylseq_task",
      perform    => 1,
      option     => "
echo 'sample,fastq_1,fastq_2,genome' > fileList1.list.csv
echo '__NAME__,__FILE__,' >> fileList1.list.csv

nextflow run $nextflow_main_nf \\
  -config $nextflow_config \\
  -profile singularity \\
  --input fileList1.list.csv \\
  --outdir . \\
  --genome $genome \\
  --igenomes_base $igenomes_base \\
  --aligner $aligner \\
  --run_targeted_sequencing \\
  --target_regions_file $convered_bed \\
  --collecthsmetrics

status=\$?
if [ \$status -ne 0 ]; then
  echo \"Error: nextflow run nf-core/methylseq failed with status \$status\"
  exit \$status
else
  echo \"nextflow run nf-core/methylseq completed successfully\"
  rm -rf work
fi
",
      program                  => "",
      check_program            => 0,
      parameterSampleFile1_ref => $source_ref,
      no_output                => 1,
      output_ext               => "bismark/summary/bismark_summary_report.html,bismark/methylation_calls/methylation_coverage/__NAME___1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
      output_to_same_folder    => 0,
      samplename_in_result     => 0,
      sh_direct                => $sh_direct,
      no_docker                => 1,
      pbs                      => {
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "50gb",
      },
    };

    push @$tasks, $nextflow_methylseq_task;

    if ($perform_downstream_analysis) {
      $methylkitprep_task = $nextflow_methylseq_task;
    }
  } ## end else [ if ( $nextflow_run_mode...)]

  if ( !$perform_downstream_analysis ) {
    return ($config);
  }

  my $methylkitcorr_task = add_MethylKitCorr( $config, $def, $tasks, $target_dir, "MethylKitCorr", [ $methylkitprep_task, ".bismark.cov.gz" ], "bismarkCoverage" );

  if ( $def->{perform_age_estimation} ) {
    my $methy_age_task = add_MethylAgeEstimation( $config, $def, $tasks, $target_dir, "dnaMethyAge", $methylkitcorr_task );
  }

  if ( $def->{pairs} ) {
    add_MethylDiffAnalysis( $config, $def, $tasks, $target_dir, $methylkitcorr_task );
  }

  $config->{"sequencetask"} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => { tasks => $tasks, },
    sh_direct  => 0,
    pbs        => {
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  return ($config);

} ## end sub getNextflowMethylationConfig


sub performNextflowMethylation {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getNextflowMethylationConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
} ## end sub performNextflowMethylation

1;
