#!/usr/bin/perl
package Pipeline::NextflowATACseq;

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

our %EXPORT_TAGS = ( 'all' => [qw(performNextflowATACseq)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';


sub initializeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType",             "FAIL" );
  initDefaultValue( $def, "cluster",               "slurm" );
  initDefaultValue( $def, "perform_preprocessing", 1 );

  return $def;
} ## end sub initializeDefaultOptions


sub getNextflowATACseqConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

  my $task_name = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $tasks = [ @$individual, @$summary ];

  my $target_dir = $def->{target_dir};

  my $igenomes_base = getValue( $def, "igenomes_base" , "" );
  my $igenomes_base_option = $igenomes_base eq "" ? "" : "--igenomes_base $igenomes_base";

  my $genome        = getValue( $def, "genome" );
  my $aligner       = getValue( $def, "aligner", "bwa" );
  my $read_length   = getValue( $def, "read_length" );

  my $nextflow_config  = getValue( $def, "nextflow_config" );
  my $nextflow_main_nf = getValue( $def, "nextflow_main_nf" );
  my $sh_direct        = getValue( $def, "sh_direct", 1 );

  my $nextflow_atacseq_option = getValue( $def, "nextflow_atacseq_option", "" );

  my $nextflow_run_mode = getValue( $def, "nextflow_run_mode", "slurm" );

  my $nextflow_atacseq_task = "nextflow_atacseq_${nextflow_run_mode}";

  $config->{$nextflow_atacseq_task} = {
    class      => "CQS::ProgramWrapper",
    target_dir => $target_dir . "/$nextflow_atacseq_task",
    option     => "
#NXF_VER=23.10.1 

rm -f __NAME__.trace.txt __NAME__.execution_report.html __NAME__.execution_timeline.html

export NXF_SYNTAX_PARSER=v1
nextflow run $nextflow_main_nf \\
  -config $nextflow_config \\
  -profile singularity \\
  --input fileList1.list.csv \\
  --outdir . \\
  --genome $genome $igenomes_base_option \\
  --aligner $aligner \\
  --read_length $read_length \\
  --narrow_peak \\
  --skip_merge_replicates \\
  --skip_consensus_peaks \\
  --skip_deseq2_qc \\
  --skip_peak_annotation \\
  --skip_plot_fingerprint \\
  --skip_plot_profile \\
  --skip_ataqv \\
  --skip_igv \\
  --max_cpus 16 --max_memory 88.GB --max_time 24.h \\
  -with-trace __NAME__.trace.txt \\
  -with-report __NAME__.execution_report.html \\
  -with-timeline __NAME__.execution_timeline.html \\
  -resume $nextflow_atacseq_option

status=\$?
if [ \$status -ne 0 ]; then
  echo \"Error: nextflow run nf-core/atacseq failed with status \$status\"
  exit \$status
else
  echo \"nextflow run nf-core/atacseq completed successfully\"
  rm -rf work
fi

# for list input files: ",
    program                             => "",
    check_program                       => 0,
    parameterSampleFile1_ref            => $source_ref,
    parameterSampleFile1_header         => "sample,fastq_1,fastq_2,replicate",
    parameterSampleFile1_join_delimiter => ",",
    parameterSampleFile1_fileFirst      => 0,
    parameterSampleFile1_col_delimiter  => ",",
    parameterSampleFile1_fileSuffix     => ".csv",
    parameterSampleFile1_suffix         => ",1",
    no_prefix                           => 1,
    no_output                           => 1,
    output_ext                          => "multiqc/bismark/multiqc_report.html",
    samplename_in_result                => 0,
    output_to_same_folder               => 1,
    sh_direct                           => $sh_direct,
    no_docker                           => 1,
    pbs                                 => {
      "nodes"    => "1:ppn=16",
      "walltime" => "24",
      "mem"      => "100gb",
    },
  };

  push @$tasks, $nextflow_atacseq_task;

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

} ## end sub getNextflowATACseqConfig


sub performNextflowATACseq {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getNextflowATACseqConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
} ## end sub performNextflowATACseq

1;
