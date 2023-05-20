#!/usr/bin/perl
package Pipeline::PIPseq;

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
use scRNA::Modules;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performPIPseq)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );
  initDefaultValue( $def, "perform_preprocessing",   0 );

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

  my $task_name = $def->{task_name};

  my $email = $def->{email};

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  my $target_dir = $def->{target_dir};

  my $pipseeker_docker_command = getValue($def, "pipseeker_docker_command");
  my $pipseeker_star_index = getValue($def, "pipseeker_star_index");
  my $singularity_option = getValue($def, "singularity_option", "");

  my $pipseeker = "pipseeker";
  $config->{$pipseeker} = {
    class => "CQS::ProgramWrapperOneToMany",
    target_dir => "${target_dir}/${pipseeker}",
    option => "

$pipseeker_docker_command count --input-path __FILE__ --id __NAME__ --output-root ${target_dir}/${pipseeker}/result/__NAME__ --star-index-path $pipseeker_star_index --star-threads 8

#__OUTPUT__
",
    interpretor => "",
    program => "",
    check_program => 0,
    source_arg => "--input-path",
    source_ref => $untrimed_ref,
    output_to_same_folder => 0,
    output_arg => "--output-root",
    samplename_in_result => 0,
    samplename_in_result => 0,
    output_file_prefix => "",
    output_file_ext => "filtered_matrix/sensitivity__ITER_",
    iteration => 5,
    iteration_fill_length => 1,
    use_tmp_folder => 0,
    sh_direct => 0,
    no_docker => 1,
    pbs => {
      "nodes" => "1:ppn=8",
      "walltime" => "96",
      "mem" => "100gb"
    },
  };

  push @$individual, ($pipseeker);

  my $pipseeker_qc = "qc";

  if(getValue($def, "perform_individual_qc", 1)){
    my $qc_pattern = getValue($def, "qc_pattern", "sensitivity_[12345]");
    add_individual_qc($config, $def, $summary, $target_dir, $pipseeker_qc, undef, [$pipseeker, $qc_pattern], undef, undef, undef);
  }

  return ($config);
}

sub performPIPseq {
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
