#!/usr/bin/perl
package MyTest::OneToOneTask;

use strict;
use warnings;
use CQS::StringUtils;
use CQS::ConfigUtils;
use CQS::Task;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_11";
  bless $self, $class;
  return $self;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );
  my $result = {};

  my $raw_files = get_raw_files($config, $section, "source");
  for my $sample_name (keys %$raw_files){
    $result->{$sample_name} = filter_array([$result_dir . "/" . $sample_name . ".csv"], $pattern);
  }
  
  return $result;
}

1;
