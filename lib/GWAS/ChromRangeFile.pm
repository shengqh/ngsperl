#!/usr/bin/perl
package GWAS::ChromRangeFile;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::UniqueTask;
use CQS::NGSCommon;
use CQS::StringUtils;
use GWAS::GwasUtils;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_cbf";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $interval_bed = get_option_file( $config, $section, "source" );
  my @chroms = readChromosomesFromBedFile( $interval_bed );
  my $lastChrom = $chroms[ scalar(@chroms) - 1 ];

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $finalFile = $task_name . "_chr" . $lastChrom . ".bed";
  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $finalFile );

  for my $chrom (@chroms) {
    print $pbs "grep \"^$chrom\\s\" > ${task_name}_chr${lastChrom}.bed \n\n";
  }
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $interval_bed = get_option_file( $config, $section, "source" );
  my @chroms = readChromosomesFromBedFile( $interval_bed );

  my $result = {};
  for my $chrom (@chroms) {
    $result->{"chr${chrom}"} = filter_array( ["$result_dir/${task_name}_chr${chrom}.bed"], $pattern );
  }
  return $result;
}

1;
