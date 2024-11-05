#!/usr/bin/perl
package BioVU::Utils;

use strict;
use warnings;
use File::Basename;
use CQS::ConfigUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(add_phenotype)] );

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
    },
    output_file_ext     => ".phenotype.html",
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
  return($config);
}

1;
