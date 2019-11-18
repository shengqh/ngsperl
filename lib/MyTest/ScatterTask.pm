#!/usr/bin/perl
package MyTest::ScatterTask;

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
use Data::Dumper;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_scatter";
  $self->{_pbskey} = "",
  bless $self, $class;
  return $self;
}

sub get_scatter_names {
  my ($self, $config, $section) = @_;
  return ["chr1", "chr2", "chr3"];  
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );
  my $result = {};

  my $scatter_names = $self->get_scatter_names($config, $section);
  for my $sample_name (@$scatter_names){
    $result->{$sample_name} = filter_array([$result_dir . "/" . $sample_name . ".csv"], $pattern);
  }
  
  return $result;
}

1;