#!/usr/bin/perl
package SmallRNA::exceRpt;

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
use CQS::ProgramWrapperOneToOne

our @ISA = qw(CQS::ProgramWrapperOneToOne);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_er";
  bless $self, $class;
  return $self;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my ( $source_files, $source_file_arg, $source_file_join_delimiter ) = get_parameter_sample_files( $config, $section, "source" );

  my $result = {};
  for my $sample_name ( sort keys %$source_files ) {
    my $sample_file = $source_files->{$sample_name}[0];
    my $sample_file_name = basename($sample_file);
    $sample_file_name =~ s/.gz$//g;

    my $cur_dir = $result_dir . "/$sample_name";

    my @result_files = "$cur_dir/${sample_file_name}_${sample_name}/EXOGENOUS_rRNA/unaligned.fq.gz";
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
