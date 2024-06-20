#!/usr/bin/perl
package Pipeline::TwistMethylation;

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

##20240617 looks like the TwistMethylation pipeline is out of date. So stop here.

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performTwistMethylation)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );
  initDefaultValue( $def, "perform_preprocessing",   1 );

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

  my $task_name = $def->{task_name};

  my $email = $def->{email};

  #$def->{perform_cutadapt} = 0;

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);
  my $tasks = [@$individual, @$summary];

  my $targetDir = $def->{target_dir};

  my $thread = getValue($def, "thread", 8);
  my $fasta_file = getValue($def, "bwameth_fasta");

  my $methylation_key = "methylation_key";

  my $bwameth_task = "bwameth";
  $config->{$bwameth_task} = {
    class      => "CQS::ProgramWrapperOneToOne",
    perform    => 1,
    program => "",
    check_program => 0,
    option     => "
bwameth.py --reference $fasta_file \\
  -t $thread \\
  --read-group '\@RG\\tSAMPLE_ID:1\\tPL:illumina\\tLB:__NAME__\\tSM:__NAME__' \\
  __FILE__ | samtools view -bS - > __NAME__.tmp.bam

status=\$?
if [ \$status -ne 0 ]; then
  echo \$status > __NAME__.bwameth.failed 
  rm -f __NAME__.bwameth.succeed __NAME__.tmp.bam
else
  rm -f __NAME__.bwameth.failed
  touch __NAME__.bwameth.succeed
  mv __NAME__.tmp.bam __NAME__.bam
fi

",
    thread     => $thread,
    target_dir => "${targetDir}/" . getNextFolderIndex($def) . "bwameth",
    parameterSampleFile1_ref => $source_ref,
    parameterSampleFile1_arg => "",
    parameterSampleFile1_join_delimiter => " ",
    output_to_same_folder => 0,
    no_output => 1,
    output_ext => ".bam",
    no_docker => 1,
    pbs        => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "24",
      "mem"       => "80gb"
    },
  };
  push(@$tasks, $bwameth_task);

  # my @report_files = ();
  # my @report_names = ();
  # my @copy_files   = ();
  # if ( defined $config->{fastqc_raw_summary} ) {
  #   push( @report_files, "fastqc_raw_summary",                   ".FastQC.baseQuality.tsv.png" );
  #   push( @report_files, "fastqc_raw_summary",                   ".FastQC.sequenceGC.tsv.png" );
  #   push( @report_files, "fastqc_raw_summary",                   ".FastQC.adapter.tsv.png" );
  #   push( @report_names, "fastqc_raw_per_base_sequence_quality", "fastqc_raw_per_sequence_gc_content", "fastqc_raw_adapter_content" );
  # }
  # if(( defined $dnmtools_task ) && ( defined $config->{$dnmtools_task} )) {
  #   push( @copy_files, $dnmtools_task, ".amr\$" );
  #   push( @copy_files, $dnmtools_task, ".hmr\$" );
  #   push( @copy_files, $dnmtools_task, ".pmd\$" );
  # }
  # if(( defined $dnmtoolsannovar_task ) && ( defined $config->{$dnmtoolsannovar_task} )) {
  #   push( @copy_files, $dnmtoolsannovar_task, ".annovar.final.tsv\$" );
  # }
  # if (( defined $methylkitcorr_task ) && ( defined $config->{$methylkitcorr_task} )) {
  #   push( @copy_files, $methylkitcorr_task, ".pdf\$|.png\$|.rds\$" );
  # }
  # if (( defined $methylkitdiff_task ) && ( defined $config->{$methylkitdiff_task} )) {
  #   push( @copy_files, $methylkitdiff_task, ".dmcpgs\$" );
  # }
  # if (( defined $methylkitdiffannovar_task ) && ( defined $config->{$methylkitdiffannovar_task} )) {
  #   push( @copy_files, $methylkitdiffannovar_task, ".annovar.final.tsv\$" );
  # }
  # if ( defined $webgestalt_task ) {
  #   push( @copy_files, $webgestalt_task, "_geneontology_Biological_Process\$" );
  #   push( @copy_files, $webgestalt_task, "_geneontology_Cellular_Component\$" );
  #   push( @copy_files, $webgestalt_task, "_geneontology_Molecular_Function\$" );
  #   push( @copy_files, $webgestalt_task, "_pathway_KEGG\$" );
  # }

  # $config->{report} = {
  #   class                      => "CQS::BuildReport",
  #   perform                    => 1,
  #   target_dir                 => "$targetDir/report",
  #   #docker_command             => "singularity exec -c -e -B /home,/gpfs51,/gpfs52,/panfs,/data,/dors,/nobackup,/tmp -H `pwd` /data/cqs/softwares/singularity/wgbs_r.1.1.sif ",
  #   docker_prefix           => "wgbs_r_",
  #   report_rmd_file            => "../Pipeline/WGBS.Rmd",
  #   additional_rmd_files       => "../Pipeline/Pipeline.R;reportFunctions.R",
  #   parameterSampleFile1_ref   => \@report_files,
  #   parameterSampleFile1_names => \@report_names,
  #   parameterSampleFile2 => {
  #     task_name => $task_name,
  #     meta_data => "../../" . $task_name . "_meta.tsv",
  #     abismal_path  => $config->{abismal}{target_dir} . "/result/",
  #     dnmtools_path => $config->{DNMTools}{target_dir} . "/result/",
  #     MethylKitCorr_path => $config->{MethylKitCorr}{target_dir} . "/result/",
  #     MethylKitDiff_path => $config->{MethylKitDiff}{target_dir} . "/result/",
  #   },
  #   parameterSampleFile3_ref   => \@copy_files,
  #   parameterSampleFile4_ref   => $webgestalt_task,
  #   parameterSampleFile5_ref   => $abismal_task,
  #   parameterSampleFile6_ref   => $dnmtools_task,
  #   parameterSampleFile7_ref   => [ $methylkitcorr_task, ".pdf\$|.png\$" ],
  #   parameterSampleFile8_ref   => $methylkitdiff_task,
  #   sh_direct                  => 1,
  #   pbs                        => {
  #     "nodes"     => "1:ppn=1",
  #     "walltime"  => "1",
  #     "mem"       => "40gb"
  #   },
  # };
  # push( @$tasks, "report" );

  $config->{"sequencetask"} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${targetDir}/sequencetask",
    option     => "",
    source     => {
      tasks => $tasks,
    },
    sh_direct => 0,
    pbs       => {
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  return($config);
};


sub performTwistMethylation {
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

1;

