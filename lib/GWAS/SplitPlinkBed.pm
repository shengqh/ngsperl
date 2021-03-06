#!/usr/bin/perl
package GWAS::SplitPlinkBed;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;
use GWAS::GwasUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_spd";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $rawFiles = get_raw_files( $config, $section );

  my $interval_bed = parse_param_file( $config, $section, "interval_bed", 1 );
  my @chroms = readChromosomesFromBedFile( $interval_bed );
  my $lastChrom = $chroms[ scalar(@chroms) - 1 ];

  for my $sampleName ( sort keys %$rawFiles ) {
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sampleName );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sampleName );
    my $log_desc = $cluster->get_log_description($log);

    my $files      = $rawFiles->{$sampleName};
    my $file       = $files->[0];
    my $filePrefix = $file;
    $filePrefix =~ s{\.[^.]+$}{};

    my $finalFile = $sampleName . ".chr" . $lastChrom . ".bed";
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $finalFile );

    for my $chrom (@chroms) {
      print $pbs "plink2 --bfile $filePrefix --chr $chrom --make-bed --out ${sampleName}_chr${chrom} \n\n";
    }
    $self->close_pbs( $pbs, $pbs_file );
  }
}

sub get_result_dependent_pbs {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $interval_bed = parse_param_file( $config, $section, "interval_bed", 1 );
  my @chroms = readChromosomesFromBedFile( $interval_bed );

  my $pbsFiles = $self->get_pbs_files($config, $section);
  
  my $result = {};
  for my $sampleName ( sort keys %$pbsFiles ) {
    my $pbs_file = $pbsFiles->{$sampleName};
    for my $chrom (@chroms) {
      $result->{$sampleName . "_chr${chrom}"} = $pbs_file;
    }
  }
  return $result;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $interval_bed = parse_param_file( $config, $section, "interval_bed", 1 );
  my @chroms = readChromosomesFromBedFile( $interval_bed );
  my $rawFiles = get_raw_files( $config, $section );

  my $result = {};
  for my $sampleName ( sort keys %$rawFiles ) {
    my @result_files = ();
    for my $chrom (@chroms) {
      $result->{$sampleName . "_chr${chrom}"} = filter_array( ["$result_dir/${sampleName}_chr${chrom}.bed"], $pattern );
    }
  }
  return $result;
}

1;
