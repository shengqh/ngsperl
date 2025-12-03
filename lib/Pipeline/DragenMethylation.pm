#!/usr/bin/perl
package Pipeline::DragenMethylation;

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

our %EXPORT_TAGS = ( 'all' => [qw(performDragenMethylation)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';


sub initializeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType",             "FAIL" );
  initDefaultValue( $def, "cluster",               "slurm" );
  initDefaultValue( $def, "perform_preprocessing", 0 );

  initDefaultValue( $def, "methylation_mincov", "20" );

  return $def;
} ## end sub initializeDefaultOptions


sub getDragenMethylationConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

  my $task_name = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $tasks = [ @$individual, @$summary ];

  my $target_dir = $def->{target_dir};

  my $multiqc_task = undef;
  if ( $def->{dragen_stats_folder} ) {
    my $dragen_stats_folder = $def->{dragen_stats_folder};
    $multiqc_task = "dragen_multiqc";
    $config->{$multiqc_task} = {
      class         => "CQS::ProgramWrapper",
      perform       => 1,
      target_dir    => "$target_dir/$multiqc_task",
      interpretor   => "",
      program       => "",
      check_program => 0,
      option        => "multiqc __FILE__
",
      parameterSampleFile1      => { $task_name => [$dragen_stats_folder] },
      parameterSampleFile1_type => "array",
      no_output                 => 1,
      samplename_in_result      => 0,
      output_file_ext           => "multiqc_report.html",
      docker_prefix             => "multiqc_",
      sh_direct                 => 0,
      pbs                       => {
        "nodes"    => "1:ppn=1",
        "walltime" => "12",
        "mem"      => "40"
      },
    };
    push( @$tasks, $multiqc_task );
  } ## end if ( $def->{dragen_stats_folder...})

  my $convert_task = "dragen_to_bismark";
  $config->{$convert_task} = {
    class       => "CQS::ProgramWrapperOneToOne",
    perform     => 1,
    target_dir  => "$target_dir/$convert_task",
    interpretor => "python",
    program     => "../Methylation/dragen_to_bismark.py",
    option      => " --input __FILE__ --output __OUTPUT__.bismark.cov.tmp.gz

status=\$?
if [ \$status -ne 0 ]; then
  echo \"dragen_to_bismark failed with status \$status\"
  rm -f __OUTPUT__.bismark.cov.tmp.gz
else
  mv __OUTPUT__.bismark.cov.tmp.gz __OUTPUT__.bismark.cov.gz
  echo \"dragen_to_bismark succeeded\"
fi
exit \$status
",
    parameterSampleFile1_ref => "files",
    output_file_ext          => ".bismark.cov.gz",
    sh_direct                => 0,
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "12",
      "mem"      => "40"
    },
  };
  push( @$tasks, $convert_task );

  my $methylkitcorr_task = add_MethylKitCorr( $config, $def, $tasks, $target_dir, "MethylKitCorr", $convert_task, "bismarkCoverage" );

  if ( $def->{perform_age_estimation} ) {
    my $methy_age_task = add_MethylAgeEstimation( $config, $def, $tasks, $target_dir, "dnaMethyAge", $methylkitcorr_task );
  }

  my $methylkitdiff_task             = undef;
  my $methylkitdiffannovar_task      = undef;
  my $MethylKitDiffAnnovarGenes_task = undef;
  my $webgestalt_task                = undef;

  if ( defined $def->{pairs} ) {
    my $task_map = add_MethylDiffAnalysis( $config, $def, $tasks, $target_dir, $methylkitcorr_task );

    #print(Dumper($task_map));

    $methylkitdiff_task             = $task_map->{methylkitdiff_task};
    $methylkitdiffannovar_task      = $task_map->{methylkitdiffannovar_task};
    $MethylKitDiffAnnovarGenes_task = $task_map->{MethylKitDiffAnnovarGenes_task};
    $webgestalt_task                = $task_map->{webgestalt_task};
  } ## end if ( defined $def->{pairs...})

  my $MethylKitCorr_path = $config->{$methylkitcorr_task}{target_dir} . "/result/";
  my $MethylKitDiff_path = ( defined $methylkitdiff_task ) ? $config->{$methylkitdiff_task}{target_dir} . "/result/" : undef;
  #print("MethylKitDiff_path=", $MethylKitDiff_path, "\n");

  my $version_files = get_version_files($config);

  my @report_files = ();
  my @report_names = ();
  my @copy_files   = ();

  if ( defined $multiqc_task ) {
    push( @report_files, $multiqc_task );
    push( @report_names, "multiqc_report_html" );
  }

  $config->{report} = {
    class      => "CQS::BuildReport",
    perform    => 1,
    target_dir => "$target_dir/report",
    #docker_command             => "singularity exec -c -e -B /home,/gpfs51,/gpfs52,/panfs,/data,/dors,/nobackup,/tmp -H `pwd` /data/cqs/softwares/singularity/wgbs_r.1.1.sif ",
    docker_prefix              => "wgbs_r_",
    report_rmd_file            => "../Pipeline/DragenMethylation.Rmd",
    additional_rmd_files       => "../Pipeline/Pipeline.R;reportFunctions.R",
    parameterSampleFile1_ref   => \@report_files,
    parameterSampleFile1_names => \@report_names,
    parameterSampleFile2       => {
      task_name            => $task_name,
      email                => $def->{email},
      affiliation          => $def->{affiliation},
      MethylKitCorr_path   => $MethylKitCorr_path,
      MethylKitDiff_path   => $MethylKitDiff_path,
      dragen_qc_file       => $def->{dragen_qc_file},
      dragen_multi_qc_html => $def->{dragen_multi_qc_html},
    },
    parameterSampleFile3_ref => \@copy_files,
    copy_to_root_folder      => 1,
    parameterSampleFile7_ref => [ $methylkitcorr_task, ".png" ],
    parameterSampleFile8_ref => [ $methylkitdiff_task, $methylkitdiffannovar_task ],
    parameterSampleFile9     => $version_files,
    sh_direct                => 1,
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "40gb"
    },
  };
  push( @$tasks, "report" );

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

} ## end sub getDragenMethylationConfig


sub performDragenMethylation {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getDragenMethylationConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
} ## end sub performDragenMethylation

1;
