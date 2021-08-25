#!/usr/bin/perl
package Annotation::FindOverlapGene;

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
  $self->{_suffix} = "_fo";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $gene_sorted_bed = get_param_file( $config->{$section}{gene_sorted_bed}, "gene_sorted_bed", 1 );

  my $script = dirname(__FILE__) . "/findOverlapGene.py";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $log = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $final_file = $self->get_final_file($config, $section, $result_dir);
  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  for my $sample_name ( sort keys %raw_files ) {
    my @bed_files = @{ $raw_files{$sample_name} };
    for my $bed_file (@bed_files) {
      my $finalFile = $bed_file . ".overlap.tsv";
      print $pbs " 
python3 $script $option -i $bed_file -o $finalFile -g $gene_sorted_bed
";
    }
  }
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my @bed_files    = @{ $raw_files{$sample_name} };
    my @result_files = ();
    for my $bed_file (@bed_files) {
      my $finalFile = $bed_file . ".overlap.tsv";

      push( @result_files, $finalFile );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
