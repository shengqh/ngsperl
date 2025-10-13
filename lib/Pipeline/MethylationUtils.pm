#!/usr/bin/perl
package Pipeline::MethylationUtils;

use strict;
use warnings;
use CQS::StringUtils;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Data::Dumper;
use List::MoreUtils qw(first_index);
use Hash::Merge     qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(
        add_MethylKitCorr
        add_MethylAgeEstimation
    )
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';


sub add_MethylKitCorr {
  my ( $config, $def, $tasks, $target_dir, $methylkitcorr_task, $methylkitprep_ref, $pipeline ) = @_;
  $config->{$methylkitcorr_task} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => "${target_dir}/" . getNextFolderIndex($def) . "MethylKitCorr",
    docker_prefix            => "wgbs_r_",
    rtemplate                => "../Methylation/methylkit_corr.R",
    rReportTemplate          => "../Methylation/methylkit_corr.Rmd;../CQS/reportFunctions.R;../CQS/countTableVisFunctions.R",
    run_rmd_independent      => 1,
    rmd_ext                  => ".methylation.corr.html",
    output_file_ext          => ".methylation.corr.html;.filtered.cpg.meth.rds",
    parameterSampleFile1_ref => $methylkitprep_ref,
    parameterSampleFile2     => {
      task_name        => getValue( $def, "task_name" ),
      email            => getValue( $def, "email" ),
      affiliation      => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
      org              => getValue( $def, "genome" ),
      var              => "group",
      control_group    => $def->{control_group},
      mincov           => getValue( $def, "methylation_mincov" ),
      corr_dim1_cutoff => $def->{corr_dim1_cutoff},
      mds_legendPos    => getValue( $def, "mds_legendPos", "bottomleft" ),
      pipeline         => $pipeline,
    },
    parameterSampleFile3_ref => "groups",
    sh_direct                => 1,
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "4",
      "mem"      => "80gb"
    },
  };
  push( @$tasks, $methylkitcorr_task );

  return ($methylkitcorr_task);
} ## end sub add_MethylKitCorr


sub add_MethylAgeEstimation {
  my ( $config, $def, $tasks, $target_dir, $methy_age_task, $methylkitcorr_task ) = @_;
  $config->{$methy_age_task} = {
    class                    => "CQS::UniqueRmd",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . "$methy_age_task",
    report_rmd_file          => "../Methylation/methylation_age.rmd",
    additional_rmd_files     => "../CQS/reportFunctions.R",
    option                   => "",
    parameterSampleFile1_ref => [ $methylkitcorr_task, 'meth.rds$' ],
    parameterSampleFile2     => {
      task_name        => getValue( $def, "task_name" ),
      email            => getValue( $def, "email" ),
      affiliation      => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
      probe_locus_file => getValue( $def, "probe_locus_file" ),
    },
    parameterSampleFile3     => $def->{"age_dict"},
    suffix                   => ".age",
    output_file_ext          => ".age.html",
    can_result_be_empty_file => 0,
    sh_direct                => 1,
    no_docker                => 1,
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "40gb"
    },
  };
  push( @$tasks, $methy_age_task );

  return ($methy_age_task);
} ## end sub add_MethylAgeEstimation

1;
