#!/usr/bin/perl
package eQTL::MatrixEQTL;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::Task;
use Pipeline::PipelineUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_me";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );
  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $snp_genotype_files    = get_raw_files( $config, $section );
  my $snp_location_files    = get_raw_files( $config, $section, "snp_pos_files" );
  my $gene_expression_files = get_raw_files( $config, $section, "gene_expression_files" );
  my $gene_location_files   = get_raw_files( $config, $section, "gene_pos_files" );
  my $prefix = get_option( $config, $section, "prefix", "" );
  my $sort_data_by_tcga = get_option( $config, $section, "sort_data_by_tcga", 0 );

  my $script = dirname(__FILE__) . "/MatrixEQTL.r";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  my $expectFiles = $self->result( $config, $section );
  for my $sampleName ( keys %$snp_genotype_files ) {
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sampleName );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sampleName );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";

    my $cis_file = "${prefix}${sampleName}.cis.txt";
    my $trans_file = "${prefix}${sampleName}.trans.txt";
    #my $final_file = basename( $expectFiles->{$sampleName}->[0] );
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $cis_file );

    my $snp_genotype_file    = $snp_genotype_files->{$sampleName}[0];
    my $snp_location_file    = $snp_location_files->{$sampleName}[0];
    my $gene_expression_file = $gene_expression_files->{$sampleName}[0];
    my $gene_location_file   = $gene_location_files->{$sampleName}[0];

    print $pbs "R --vanilla -f $script --args $snp_genotype_file $snp_location_file $gene_expression_file $gene_location_file $sort_data_by_tcga $cis_file $trans_file\n";
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section, 0 );

  my $snp_genotype_files = get_raw_files( $config, $section );
  my $prefix = get_option( $config, $section, "prefix", "" );

  my $result = {};

  for my $sampleName (keys %$snp_genotype_files) {
    my @result_files = ();
    push( @result_files, "${result_dir}/${prefix}${sampleName}.cis.txt" );
    push( @result_files, "${result_dir}/${prefix}${sampleName}.trans.txt" );
    $result->{$sampleName} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
