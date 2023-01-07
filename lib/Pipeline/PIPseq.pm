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

  my $pipseeker_sif = getValue($def, "pipseeker_sif");
  my $pipseeker_star_index = getValue($def, "pipseeker_star_index");
  my $singularity_option = getValue($def, "singularity_option", "");

  my $pipseeker = "pipseeker";
  $config->{$pipseeker} = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "${target_dir}/${pipseeker}",
    option                => "

singularity run $singularity_option $pipseeker_sif count --input-path __FILE__ --id __NAME__ --output-root ${target_dir}/${pipseeker}/result/__NAME__ --star-index-path $pipseeker_star_index --star-threads 8

#__OUTPUT__
",
    interpretor           => "",
    program               => "",
    check_program         => 0,
    source_arg            => "-i",
    source_ref            => $untrimed_ref,
    output_to_same_folder => 0,
    output_arg            => "--output-root",
    output_no_name => 1,
    output_file_ext    => "filtered_matrix/sensitivity_1",
    use_tmp_folder        => 0,
    sh_direct             => 0,
    pbs                   => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "96",
      "mem"       => "100gb"
    },
  };
  #push @$individual, ($cutruntools2);

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
