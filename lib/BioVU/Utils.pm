#!/usr/bin/perl
package BioVU::Utils;

use strict;
use warnings;
use File::Basename;
use CQS::ConfigUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(add_phenotype add_linear_association)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.02';

sub add_phenotype {
  my ($config, $def, $task_name, $phecode_dic) = @_;

  my $report_file = dirname(__FILE__) . "/../CQS/reportFunctions.R";
  my $rmd_file = dirname(__FILE__) . "/prepare_phenotype_data.rmd";

  $config->{$task_name} = {
    class                => "CQS::ProgramWrapperOneToOne",
    perform              => 1,
    target_dir           => $def->{target_dir} . "/" . $task_name,
    program              => "",
    check_program        => 0,
    option               => "
cp -f $report_file .    
cp -f $rmd_file __NAME__.phenotype.rmd
echo -e '__NAME__\\tphename\\n__FILE__\\tphecode\\n../__FILE2__\\toption_file' > input_options.txt
R -e \"library(knitr);rmarkdown::render(input='__NAME__.phenotype.rmd');\"   
",
    parameterSampleFile1 => $phecode_dic,
    parameterSampleFile2 => {
      agd_file => getValue($def, "agd_file"),
      phecode_data_file => getValue($def, "phecode_data_file"),
      phecode_map_file => getValue($def, "phecode_map_file"),
      min_occurance => getValue($def, "min_occurance", 1),
      email => getValue($def, "email", ""),
      affiliation => getValue($def, "affiliation", ""),
      show_code => getValue($def, "show_code", 1),
    },
    output_file_ext     => ".phenotype.csv",
    output_to_same_folder => 0,
    sh_direct            => 1,
    no_prefix => 1,
    no_output           => 1,
    pbs                  => {
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "20gb"
    },
  };

  my $summary_task = $task_name . "_summary";
  $config->{$summary_task} = {
    class => "CQS::UniqueR",
    perform => 1,
    target_dir => $def->{target_dir} . "/" . $summary_task,
    option => "",
    rtemplate => "../BioVU/prepare_phenotype_data_summary.r",
    parameterSampleFile1_ref => $task_name,
    output_file_ext => ".phenotype_summary.csv",
    sh_direct => 1,
    no_output => 1,
    pbs => {
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "5gb"
    },
  };

  return($config);
}

sub add_linear_association {
  my ($config, $def, $task_name, $phecode_dic, $phenotype_ref) = @_;

  my $report_file = dirname(__FILE__) . "/../CQS/reportFunctions.R";
  my $rmd_file = dirname(__FILE__) . "/linear_association.rmd";

  $config->{$task_name} = {
    class                => "CQS::ProgramWrapperOneToOne",
    perform              => 1,
    target_dir           => $def->{target_dir} . "/" . $task_name,
    program              => "",
    check_program        => 0,
    option               => "
cp -f $report_file .    

cp -f $rmd_file __NAME__.linear_association.rmd

echo -e '__NAME__\\tphename\\n__FILE__\\tphecode\\n__FILE3__\\tphefile' > input_options.txt

R -e \"library(knitr);rmarkdown::render(input='__NAME__.linear_association.rmd');\"   

#__FILE2__
",
    parameterSampleFile1 => $phecode_dic,
    parameterSampleFile2 => {
      genotype_file => getValue($def, "genotype_file"),
      pca_file => getValue($def, "pca_file"),
      ancestry_file => getValue($def, "ancestry_file"),
      age_file => getValue($def, "age_file"),
      demographics_file => getValue($def, "demographics_file"),
      phecode_map_file => getValue($def, "phecode_map_file"),
      email => getValue($def, "email", ""),
      affiliation => getValue($def, "affiliation", ""),
    },
    parameterSampleFile3_ref => $phenotype_ref,
    output_to_same_folder => 0,
    output_file_ext => ".linear_association.csv",
    no_prefix => 1,
    sh_direct            => 0,
    no_output           => 1,
    pbs                  => {
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "20gb"
    },
  };

  my $summary_task = $task_name . "_summary";
  $config->{$summary_task} = {
    class => "CQS::UniqueR",
    perform => 1,
    target_dir => $def->{target_dir} . "/" . $summary_task,
    option => "",
    rtemplate => "../BioVU/linear_association_summary.r",
    parameterSampleFile1_ref => $task_name,
    parameterFile1 => getValue($def, "phecode_map_file"),
    output_file_ext => ".linear_association.csv",
    sh_direct => 0,
    no_output => 1,
    pbs => {
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "5gb"
    },
  },
};

1;
