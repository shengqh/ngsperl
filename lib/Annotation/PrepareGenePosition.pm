#!/usr/bin/perl
package Annotation::PrepareGenePosition;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_pg";
  $self->{_docker_prefix} = "genepos_";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $datasetName = get_option( $config, $section, "dataset_name" );
  my $geneNames   = get_option( $config, $section, "gene_names" );
  my $isGFF       = get_option( $config, $section, "output_gff" );
  my $addChr       = get_option( $config, $section, "add_chr" );

  my $script = dirname(__FILE__) . "/prepareGenePosition.r";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $log = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $finalFile = $isGFF ? $task_name . ".genes.gff" : $task_name . ".genes.tsv";
  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $finalFile );

  print $pbs " 
R --vanilla -f $script --args $datasetName $finalFile $isGFF $addChr $geneNames
";
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $isGFF = get_option( $config, $section, "output_gff" );
  my $finalFile = $isGFF ? $task_name . ".genes.gff" : $task_name . ".genes.tsv";

  my $result       = {};
  my @result_files = ();
  push( @result_files, "$result_dir/$finalFile" );
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
