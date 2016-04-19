#!/usr/bin/perl
package CQS::UniqueTask;

use strict;
use warnings;
use CQS::Task;
use CQS::ConfigUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_pbskey} = "";
  bless $self, $class;
  return $self;
}

sub get_clear_map {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $comparison_result = $self->SUPER::get_clear_map( $config, $section, $pattern );
  my $result            = {};
  my @result_files      = ();
  for my $crs ( sort values %$comparison_result ) {
    for my $cr (@$crs) {
      push( @result_files, $cr );
    }
  }
  $result->{$task_name} = \@result_files;

  return $result;
}

1;
