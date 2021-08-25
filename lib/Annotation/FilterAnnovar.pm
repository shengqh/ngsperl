#!/usr/bin/perl
package Annotation::FilterAnnovar;

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
  $self->{_suffix} = "_fa";
  bless $self, $class;
  return $self;
}

sub get_maximum_freq_values {
  my ( $config, $section ) = @_;
  my $exac_values = get_option( $config, $section, "maximum_freq_values", "0.01,0.001" );
  my @result = split( ',', $exac_values );
  return @result;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my @freq_values = get_maximum_freq_values( $config, $section );
  my $sampleNamePattern = get_option( $config, $section, "sample_name_pattern", "" );
  my $sampleNameSuffix = "";
  if ( $sampleNamePattern ne "" ) {
    $sampleNameSuffix = get_option( $config, $section, "sample_name_suffix" );
    $sampleNamePattern = "-r $sampleNamePattern";
  }

  my $exac_key = get_option( $config, $section, "exac_key", "" );
  if ( $exac_key ne "" ) {
    $exac_key = "--exac_key $exac_key";
  }

  my $g1000_key = get_option( $config, $section, "g1000_key", "" );
  if ( $g1000_key ne "" ) {
    $g1000_key = "--g1000_key $g1000_key";
  }

  my $gnomad_key = get_option( $config, $section, "gnomad_key", "" );
  if ( $gnomad_key ne "" ) {
    $gnomad_key = "--gnomad_key $gnomad_key";
  }

  my $filter_fq_equal_1 = get_option( $config, $section, "filter_fq_equal_1", 0 );
  if ( $filter_fq_equal_1 ) {
    $filter_fq_equal_1 = "--filter_fq_equal_1";
  }else{
    $filter_fq_equal_1 = "";
  }

  my $script = dirname(__FILE__) . "/filterAnnovar.py";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $log = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $final_file = $self->get_final_file($config, $section, $result_dir);
  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  for my $sample_name ( sort keys %raw_files ) {
    my @annovar_files = @{ $raw_files{$sample_name} };
    my $annovar_file  = $annovar_files[0];

    if ( scalar(@freq_values) > 0 ) {
      for my $freq_value (@freq_values) {
        my $finalFilePrefix = "${sample_name}${sampleNameSuffix}.freq${freq_value}";
        my $finalFile       = $finalFilePrefix . ".gene.missense.tsv";
        print $pbs " 
python3 $script $option -i $annovar_file -t $freq_value -o $finalFilePrefix $exac_key $g1000_key $gnomad_key $filter_fq_equal_1 $sampleNamePattern 
";
      }
    }
    else {
      my $finalFilePrefix = "${sample_name}${sampleNameSuffix}";
      my $finalFile       = $finalFilePrefix . ".gene.missense.tsv";
      print $pbs " 
python3 $script $option -i $annovar_file -o $finalFilePrefix $sampleNamePattern
";
    }
  }
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my @freq_values = get_maximum_freq_values( $config, $section );
  my $sampleNamePattern = get_option( $config, $section, "sample_name_pattern", "" );
  my $sampleNameSuffix = "";
  if ( $sampleNamePattern ne "" ) {
    $sampleNameSuffix = get_option( $config, $section, "sample_name_suffix" );
  }

  my $result = {};

  for my $sample_name ( sort keys %raw_files ) {
    my @result_files = ();
    if ( scalar(@freq_values) > 0 ) {
      for my $freq_value (@freq_values) {
        push( @result_files, "$result_dir/${task_name}${sampleNameSuffix}.freq${freq_value}.filtered.tsv" );
        push( @result_files, "$result_dir/${task_name}${sampleNameSuffix}.freq${freq_value}.filtered.missense.tsv" );
        push( @result_files, "$result_dir/${task_name}${sampleNameSuffix}.freq${freq_value}.snv.missense.tsv" );
        push( @result_files, "$result_dir/${task_name}${sampleNameSuffix}.freq${freq_value}.gene.missense.tsv" );
      }
    }
    else {
      push( @result_files, "$result_dir/${task_name}${sampleNameSuffix}.filtered.tsv" );
      push( @result_files, "$result_dir/${task_name}${sampleNameSuffix}.filtered.missense.tsv" );
      push( @result_files, "$result_dir/${task_name}${sampleNameSuffix}.snv.missense.tsv" );
      push( @result_files, "$result_dir/${task_name}${sampleNameSuffix}.gene.missense.tsv" );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
